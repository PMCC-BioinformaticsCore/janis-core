

from abc import ABC, abstractmethod
from typing import Any, Optional

from janis_core.ingestion.galaxy import expressions
from janis_core.ingestion.galaxy.expressions.patterns import COMPOUND_OPT

from ...model import XMLParam
from ..tokens import Token
from ..tokens import TokenType
from ..components import Tee
from ..components import StreamMerge
from ..components import factory
from .ExecutionPath import EPathPosition
from . import analysis


class Annotator(ABC):

    """extracts a command component at the current position in an epath"""
    def __init__(self, ptr: int, positions: list[EPathPosition]):
        self.ptr = ptr
        self.positions = positions

    def update_epath_components(self, ptr: int, component: Any) -> None:
        self.positions[ptr].component = component

    @abstractmethod
    def can_annotate(self) -> bool:
        """confirm whether this type of component should be extracted at the current position"""
        ...
    
    @abstractmethod
    def do_annotate(self) -> None:
        """perform the necessary operations now that the component is valid"""
        ...
   
    @abstractmethod
    def calculate_next_ptr_pos(self) -> int:
        """update the ptr pos to continue iteration. based on component span"""
        ...
    

# --------------
# TOOL ARGUMENTS
# --------------

class OptionAnnotator(Annotator):
    """
    this one is complex, others are straightforward.
    identifies options with the following forms:
        - prefix value
        - prefix value value...
        - prefix sep value
        - prefix sep value value...

    """
    def __init__(self, ptr: int, positions: list[EPathPosition]):
        super().__init__(ptr, positions)
        self.stop_ptr: int = 0

    def can_annotate(self) -> bool:
        cpos = self.positions[self.ptr]
        npos = self.positions[self.ptr + 1]
        if analysis.is_option(cpos, npos):
            if not analysis.has_compound_structure(cpos, npos):
                return True
        return False
    
    def do_annotate(self) -> None:
        cpos = self.positions[self.ptr]
        npos = self.positions[self.ptr + 1]
        
        values_span = analysis.get_option_values_span(self.ptr, self.positions)
        component_span = [self.ptr, values_span[1]]

        if analysis.has_kvlinker_structure(cpos, npos):
            prefix = cpos.token.text
            separator = npos.token.text
            values = [self.positions[self.ptr + 2].token.text]

        elif analysis.has_prefix_value_structure(cpos, npos):
            prefix = cpos.token.text
            separator = ' '
            values = [pos.token.text for pos in self.positions[values_span[0]:values_span[1] + 1]]

        gxparams = [pos.token.gxparam for pos in self.positions[component_span[0]:component_span[1] + 1] if pos.token.gxparam]
        unique_params = set([x.name for x in gxparams])
        if len(unique_params) > 1:
            raise RuntimeError
        if len(gxparams) == 0:
            param = None
        elif len(gxparams) >= 1:
            param = gxparams[0]
 
        option = factory.option(
            prefix=prefix,
            separator=separator,
            values=values,
            gxparam=param
        )

        for ptr in range(component_span[0], component_span[1] + 1):
            self.update_epath_components(ptr, option)

        self.ptr = component_span[1]

    def calculate_next_ptr_pos(self) -> int:
        return self.ptr + 1



class CompoundOptionAnnotator(Annotator):

    def can_annotate(self) -> bool:
        cpos = self.positions[self.ptr]
        npos = self.positions[self.ptr + 1]
        if analysis.is_option(cpos, npos):
            if analysis.has_compound_structure(cpos, npos):
                return True
        return False
    
    def do_annotate(self) -> None:
        cpos = self.positions[self.ptr] 
        match = expressions.get_matches(cpos.token.text, COMPOUND_OPT)[0]
        option = factory.option(match.group(1), gxparam=cpos.token.gxparam, separator='', values=[match.group(2)])
        self.update_epath_components(self.ptr, option)
    
    def calculate_next_ptr_pos(self) -> int:
        return self.ptr + 1



class FlagAnnotator(Annotator):

    def can_annotate(self) -> bool:
        cpos = self.positions[self.ptr]
        npos = self.positions[self.ptr + 1]
        if analysis.is_flag(cpos, npos):
            return True
        return False
    
    def do_annotate(self) -> None:
        cpos = self.positions[self.ptr]
        flag = factory.flag(prefix=cpos.token.text, gxparam=cpos.token.gxparam)
        self.update_epath_components(self.ptr, flag)

    def calculate_next_ptr_pos(self) -> int:
        return self.ptr + 1


class PositionalAnnotator(Annotator):

    BANNED_POSITIONALS_TEXT = set([
        '', ' ', 'true', 'false', '=', ':'
    ])

    def can_annotate(self) -> bool:
        cpos = self.positions[self.ptr]
        if not analysis.is_positional(cpos):
            return False
        if cpos.token.text in self.BANNED_POSITIONALS_TEXT:
            return False
        return True
    
    def do_annotate(self) -> None:
        cpos = self.positions[self.ptr]
        positional = factory.positional(cpos.token.text, cpos.token.gxparam)
        self.update_epath_components(self.ptr, positional)
    
    def calculate_next_ptr_pos(self) -> int:
        return self.ptr + 1
 

# --------------
# LINUX CONSTRUCTS
# --------------

class StreamMergeAnnotator(Annotator):

    def can_annotate(self) -> bool:
        cpos = self.positions[self.ptr]
        if cpos.token.ttype == TokenType.LINUX_STREAM_MERGE:
            return True
        return False
    
    def do_annotate(self) -> None:
        cpos = self.positions[self.ptr]
        component = StreamMerge(cpos.token.text)
        self.update_epath_components(self.ptr, component)
   
    def calculate_next_ptr_pos(self) -> int:
        return self.ptr + 1


class RedirectAnnotator(Annotator):

    def can_annotate(self) -> bool:
        cpos = self.positions[self.ptr]
        npos = self.positions[self.ptr + 1]
        if cpos.token.ttype == TokenType.LINUX_REDIRECT and not npos.token.ttype == TokenType.END_STATEMENT:
            return True
        return False
    
    def do_annotate(self) -> None:
        cpos = self.positions[self.ptr]
        npos = self.positions[self.ptr + 1]
        component = factory.redirect_output(cpos.token.text, npos.token.text, gxparam=npos.token.gxparam)
        self.update_epath_components(self.ptr, component)
        self.update_epath_components(self.ptr + 1, component)
   
    def calculate_next_ptr_pos(self) -> int:
        return self.ptr + 2


class TeeAnnotator(Annotator):

    def can_annotate(self) -> bool:
        cpos = self.positions[self.ptr]
        npos = self.positions[self.ptr + 1]
        if cpos.token.ttype == TokenType.LINUX_TEE and not npos.token.ttype == TokenType.END_STATEMENT:
            return True
        return False

    def do_annotate(self) -> None:
        tee = Tee()
        ptr = self.ptr + 1 

        # tee arguments
        token = self.positions[ptr].token
        while token.text.startswith('-'):
            tee.options.append(token.text)
            self.update_epath_components(ptr, tee)
            ptr += 1
            token = self.positions[ptr].token

        # tee files: consumes to end of statement
        while ptr < len(self.positions) - 1:
            token = self.positions[ptr].token
            tee.files.append(token.text)
            self.update_epath_components(ptr, tee)
            ptr += 1

    def calculate_next_ptr_pos(self) -> int:
        return len(self.positions)


