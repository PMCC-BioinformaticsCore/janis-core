

from typing import Type

from janis_core.ingestion.galaxy.gxtool.model import XMLTool
from ..tokens import TokenType
from ..Command import Command
from .ExecutionPath import ExecutionPath
from .Annotator import (
    Annotator, 
    CompoundOptionAnnotator, 
    FlagAnnotator, 
    OptionAnnotator, 
    PositionalAnnotator, 
    RedirectAnnotator, 
    StreamMergeAnnotator, 
    TeeAnnotator
)

from . import analysis


"""
iterates through a ExecutionPath, yielding the current tokens being assessed.
keeps track of the location we are in the ExecutionPath

Example ExecutionPath:
abricate $sample_name --no-header --minid=80 --db = card

0            1             2            3           4          5        6           7          8
abricate     $sample_name  --no-header  --minid     =          80       --db        =          card
RAW_STRING   ENV_VAR       RAW_STRING   RAW_STRING  KV_LINKER  RAW_INT  RAW_STRING  KV_LINKER  RAW_STRING
positional1  positional2   flag1        option1     option1    option1  option2     option2    option2

at token 2, we would return {
    ctoken: token 2  (--no-header's token),
    ntokens: []
}
at token 3, we would return {
    ctoken: token 3 (--minid's token),
    ntokens: [token 4, token 5]
}

"""

linux_constructs: list[Type[Annotator]] = [
    StreamMergeAnnotator,  
    RedirectAnnotator,  
    TeeAnnotator,  
]

# annotation order matters 
tool_arguments: list[Type[Annotator]] = [
    OptionAnnotator,            # priority 1
    CompoundOptionAnnotator,    # priority 2
    FlagAnnotator,              # priority 3
    PositionalAnnotator         # priority 4
]       


class GreedyExecutionPathAnnotator:
    def __init__(self, epath: ExecutionPath, xmltool: XMLTool, command: Command):
        self.epath = epath 
        self.xmltool = xmltool
        self.command = command

    def annotate(self) -> ExecutionPath:
        self.annotate_existing_components()
        self.annotate_linux_constructs()
        self.annotate_tool_arguments()
        return self.epath

    def annotate_existing_components(self) -> None:
        self.annotate_existing_flags()
        self.annotate_existing_options()

    def annotate_existing_flags(self) -> None:
        for pos in self.epath.positions:
            if pos.token.text in self.command.flags and not pos.token.text in self.command.options:
                flag = self.command.flags[pos.token.text] 
                if flag.confidence.value == 3:
                    pos.component = self.command.flags[pos.token.text]

    def annotate_existing_options(self) -> None:
        # TODO TESTS
        ptr = 0
        while ptr < len(self.epath.positions) - 1:
            pos = self.epath.positions[ptr]

            # do we already have an option with this prefix?
            if pos.token.text in self.command.options and not pos.token.text in self.command.flags:
                option = self.command.options[pos.token.text]
                
                # if option with high confidence, annotate epath positions with this option
                if option.confidence.value == 3:
                    start, stop = analysis.get_option_values_span(ptr, self.epath.positions)
                    for i in range(ptr, stop + 1):
                        pos = self.epath.positions[i]
                        pos.component = option
                    ptr = stop
                                        
            ptr += 1

    def annotate_linux_constructs(self) -> None:
        self.iter_annotate(annotators=linux_constructs)
    
    def annotate_tool_arguments(self) -> None:
        self.iter_annotate(annotators=tool_arguments)

    def iter_annotate(self, annotators: list[Type[Annotator]]) -> None:
        ptr = 0 # reset
        while ptr < len(self.epath.positions) - 1:
            if self.can_annotate_position(ptr):
                ptr = self.do_annotate_position(ptr, annotator_classes=annotators)
            else:
                ptr += 1

    def can_annotate_position(self, ptr: int) -> bool:
        ignore_tokens = [
            TokenType.END_STATEMENT,
            TokenType.FUNCTION_CALL,
            TokenType.BACKTICK_SHELL_STATEMENT
        ]
        pos = self.epath.positions[ptr]
        if pos.token.ttype in ignore_tokens:
            return False
        if pos.component:
            return False
        return True

    def do_annotate_position(self, ptr: int, annotator_classes: list[Type[Annotator]]) -> int:
        for annotator_class in annotator_classes:
            ann = annotator_class(ptr, self.epath.positions)
            if ann.can_annotate():
                ann.do_annotate()
                return ann.calculate_next_ptr_pos()
        return ptr + 1



