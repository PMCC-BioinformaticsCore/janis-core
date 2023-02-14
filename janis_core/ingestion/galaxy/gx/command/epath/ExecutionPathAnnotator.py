

from typing import Type

from ...gxtool.tool import XMLToolDefinition
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

# annotation order matters ?
tool_arguments: list[Type[Annotator]] = [
    OptionAnnotator,   # priority 1
    CompoundOptionAnnotator,
    FlagAnnotator,
    PositionalAnnotator
]


class GreedyExecutionPathAnnotator:
    def __init__(self, epath: ExecutionPath, xmltool: XMLToolDefinition, command: Command):
        self.epath = epath 
        self.xmltool = xmltool
        self.command = command

    def annotate(self) -> ExecutionPath:
        self.mark_ignore_tokens()
        self.annotate_via_param_args()
        self.annotate_linux_constructs()
        self.annotate_tool_arguments()
        self.transfer_gxparams()
        return self.epath

    def mark_ignore_tokens(self) -> None:
        ignore_tokens = [
            TokenType.FUNCTION_CALL,
            TokenType.BACKTICK_SHELL_STATEMENT,
        ]
        for position in self.epath.positions:
            if position.token.ttype in ignore_tokens:
                position.ignore = True
        
    def annotate_via_param_args(self) -> None:
        arguments: set[str] = set([param.argument for param in self.xmltool.inputs.list() if param.argument]) # type: ignore
        ptr = 0
        while ptr < len(self.epath.positions) - 1:
            token = self.epath.positions[ptr].token
            if token.text in arguments:
                token.ttype = TokenType.FORCED_PREFIX
                ptr = self.annotate_position(ptr, annotators=tool_arguments)
            else:
                ptr += 1
    
    def annotate_linux_constructs(self) -> None:
        self.iter_annotate(annotators=linux_constructs)
    
    def annotate_tool_arguments(self) -> None:
        self.iter_annotate(annotators=tool_arguments)

    def iter_annotate(self, annotators: list[Type[Annotator]]) -> None:
        ptr = 0 # reset
        while ptr < len(self.epath.positions) - 1:
            position = self.epath.positions[ptr]
            if not position.ignore:
                ptr = self.annotate_position(ptr, annotators=annotators)
            else:
                ptr += 1

    def annotate_position(self, ptr: int, annotators: list[Type[Annotator]]) -> int:
        for annotator in annotators:
            a = annotator(ptr, self.epath.positions)
            a.annotate()
            if a.success:
                return a.calculate_next_ptr_pos()
        return ptr + 1

    def transfer_gxparams(self) -> None:
        for position in self.epath.positions:
            if position.component and position.token.gxparam:
                position.component.gxparam = position.token.gxparam

