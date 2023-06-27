

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
from . import utils as component_utils


class Annotator(ABC):

    """extracts a command component at the current position in an epath"""
    def __init__(self, ptr: int, positions: list[EPathPosition]):
        self.ptr = ptr
        self.positions = positions

        self.ctoken: Token = self.positions[self.ptr].token
        self.ntoken: Token = self.positions[self.ptr + 1].token
        self.success: bool = False

    def annotate(self) -> None:
        if self.passes_check():
            self.handle()
            self.success = True

    def update_epath_components(self, ptr: int, component: Any) -> None:
        self.positions[ptr].component = component
        self.positions[ptr].ignore = True

    @abstractmethod
    def passes_check(self) -> bool:
        """confirm whether this type of component should be extracted at the current position"""
        ...
    
    @abstractmethod
    def handle(self) -> None:
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

    """this one is complex, others are straightforward"""
    def __init__(self, ptr: int, positions: list[EPathPosition]):
        super().__init__(ptr, positions)
        self.stop_ptr: int = 0

    def passes_check(self) -> bool:
        # if self.ctoken.text == '--make_plots':
        #     print()
        if component_utils.is_option(self.ctoken, self.ntoken):
            if not component_utils.has_compound_structure(self.ctoken):
                return True
        return False
    
    def handle(self) -> None:
        value_tokens = self.get_option_values()
        gxparam = self.get_param(value_tokens)
        option = factory.option(
            prefix=self.ctoken.text,
            gxparam=gxparam,
            separator=self.get_delim(),
            values=[v.text for v in value_tokens]
        )
        for pos in range(self.ptr, self.stop_ptr + 1):
            self.update_epath_components(pos, option)

    def get_param(self, value_tokens: list[Token]) -> Optional[XMLParam]:
        out: Optional[XMLParam] = self.ctoken.gxparam
        if not out:
            for token in value_tokens:
                if token.gxparam:
                    out = token.gxparam
                    break
        return out

    def get_delim(self) -> str:
        if self.ntoken.ttype == TokenType.KV_LINKER:
            return self.ntoken.text
        return ' ' # fallback default

    def get_option_values(self) -> list[Token]:
        if self.is_kv_pair():
            return self.get_single_value()
        else:
            return self.get_multiple_values()

    def is_kv_pair(self) -> bool:
        if self.positions[self.ptr + 1].token.ttype == TokenType.KV_LINKER:
            return True
        return False
     
    def get_single_value(self) -> list[Token]:
        # if is_kv_pair() -- self.ptr + 1 will always be KV_LINKER
        self.stop_ptr = self.ptr + 2
        value = self.positions[self.ptr + 2]
        return [value.token]

    def get_multiple_values(self) -> list[Token]:
        out: list[Token] = []

        # define a max possible end point to consume upto (handles edge cases)
        final_pos = self.get_greedy_consumption_end()  
        values_type = self.positions[self.ptr + 1].token.ttype

        # look at next ntoken to see its probably a value for the option
        curr_pos = self.ptr + 1
        while curr_pos <= final_pos:
            if self.should_eat_value(curr_pos, values_type):
                out.append(self.positions[curr_pos].token)
                curr_pos += 1 
            else:
                break
        self.stop_ptr = curr_pos - 1
        return out

    def should_eat_value(self, curr_pos: int, values_type: TokenType) -> bool:
        ctoken = self.positions[curr_pos].token
        ntoken = self.positions[curr_pos + 1].token
        if component_utils.is_positional(ctoken) and ctoken.ttype == values_type and ntoken.ttype != TokenType.KV_LINKER:
            return True
        return False

    def get_greedy_consumption_end(self) -> int:
        """
        define an end point -> tee or redirect or end of positions
        most cases: end point = len(positions) to the end, first tee, first redirect etc - 2
        special case: 2nd last position. --files $input1 set end point = the above - 1
        """
        final_pos = self.get_end_token_position()
        if self.ptr < final_pos - 1:
            final_pos -= 1
        return final_pos

    def get_end_token_position(self) -> int:
        end_tokens = [TokenType.END_STATEMENT, TokenType.LINUX_TEE, TokenType.LINUX_REDIRECT]
        for position in self.positions:
            if position.token.ttype in end_tokens:
                return position.ptr - 1
        return len(self.positions) - 2  # this will never happen
    
    def calculate_next_ptr_pos(self) -> int:
        return self.stop_ptr + 1


class CompoundOptionAnnotator(Annotator):

    def passes_check(self) -> bool:
        if component_utils.is_option(self.ctoken, self.ntoken):
            if component_utils.has_compound_structure(self.ctoken):
                return True
        return False
    
    def handle(self) -> None:
        match = expressions.get_matches(self.ctoken.text, COMPOUND_OPT)[0]
        option = factory.option(match.group(1), gxparam=self.ctoken.gxparam, separator='', values=[match.group(2)])
        self.update_epath_components(self.ptr, option)
    
    def calculate_next_ptr_pos(self) -> int:
        return self.ptr + 1


class FlagAnnotator(Annotator):

    def passes_check(self) -> bool:
        if component_utils.is_flag(self.ctoken, self.ntoken):
            return True
        return False
    
    def handle(self) -> None:
        prefix = self.ctoken.text
        gxparam = self.ctoken.gxparam
        flag = factory.flag(prefix, gxparam)
        self.update_epath_components(self.ptr, flag)

    def calculate_next_ptr_pos(self) -> int:
        return self.ptr + 1


class PositionalAnnotator(Annotator):

    def passes_check(self) -> bool:
        if component_utils.is_positional(self.ctoken):
            return True
        return False
    
    def handle(self) -> None:
        value = self.ctoken.text
        gxparam = self.ctoken.gxparam
        positional = factory.positional(value, gxparam)
        self.update_epath_components(self.ptr, positional)
    
    def calculate_next_ptr_pos(self) -> int:
        return self.ptr + 1
 

# --------------
# LINUX CONSTRUCTS
# --------------

class StreamMergeAnnotator(Annotator):

    def passes_check(self) -> bool:
        if self.ctoken.ttype == TokenType.LINUX_STREAM_MERGE:
            return True
        return False
    
    def handle(self) -> None:
        component = StreamMerge(self.ctoken.text)
        self.update_epath_components(self.ptr, component)
   
    def calculate_next_ptr_pos(self) -> int:
        return self.ptr + 1


class RedirectAnnotator(Annotator):

    def passes_check(self) -> bool:
        if self.ctoken.ttype == TokenType.LINUX_REDIRECT:
            return True
        return False
    
    def handle(self) -> None:
        if self.ctoken and self.ntoken:
            component = factory.redirect_output(self.ctoken.text, self.ntoken.text, gxparam=self.ntoken.gxparam)
            self.update_epath_components(self.ptr, component)
            self.update_epath_components(self.ptr + 1, component)
   
    def calculate_next_ptr_pos(self) -> int:
        return self.ptr + 2


class TeeAnnotator(Annotator):

    def passes_check(self) -> bool:
        if self.ctoken.ttype == TokenType.LINUX_TEE:
            return True
        return False

    def handle(self) -> None:
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


