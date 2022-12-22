

from typing import Any
from enum import Enum, auto

from janis_core.types import File, DataType, Stdout
from janis_core.utils.secondary import apply_secondary_file_format_to_filename
from janis_core import (
    ToolOutput, 
    CommandTool,
    WildcardSelector,
    InputSelector,
    Selector,
    Filename
)

from ...entity_trace import trace_janis_entities
from ...unwrap import unwrap_expression
from ... import nfgen_utils
from ... import secondaries
from .. import inputs

from .model import (
    ProcessOutput,
    PathProcessOutput,
    ValProcessOutput,
    TupleProcessOutput,
    StdoutProcessOutput,
    SecondariesArrayProcessOutput
)


class OType(Enum):
    STDOUT              = auto()
    NON_FILE            = auto()
    FILE                = auto()
    FILE_PAIR           = auto()
    FILE_ARRAY          = auto()
    SECONDARIES         = auto()
    SECONDARIES_ARRAY   = auto()


def get_otype(out: ToolOutput) -> OType:
    if is_stdout_type(out):
        return OType.STDOUT
    
    elif is_secondary_type(out) and is_array_type(out):
        return OType.SECONDARIES_ARRAY
    
    elif is_secondary_type(out):
        return OType.SECONDARIES
    
    elif is_file_type(out) and is_array_type(out) and has_n_collectors(out, n=2):
        return OType.FILE_PAIR
    
    elif is_file_type(out) and is_array_type(out) and has_n_collectors(out, n=1):
        return OType.FILE_ARRAY
    
    elif is_file_type(out):
        return OType.FILE
    
    elif is_non_file_type(out) and has_n_collectors(out, n=1):
        return OType.NON_FILE

    else:
        raise NotImplementedError


def is_stdout_type(out: ToolOutput) -> bool:
    if isinstance(out.output_type, Stdout):
        return True
    return False

def is_file_type(out: ToolOutput) -> bool:
    basetype = nfgen_utils.get_base_type(out.output_type)
    if isinstance(basetype, File):
        return True
    return False

def is_secondary_type(out: ToolOutput) -> bool:
    if secondaries.is_secondary_type(out.output_type):
        return True
    return False

def is_array_type(out: ToolOutput) -> bool:
    if out.output_type.is_array():
        return True
    return False

def is_non_file_type(out: ToolOutput) -> bool:
    basetype = nfgen_utils.get_base_type(out.output_type)
    if not isinstance(basetype, File):
        return True
    return False

def has_n_collectors(out: ToolOutput, n: int) -> bool:
    sel = out.selector
    if sel is None:
        collectors = 0
    elif isinstance(sel, str):
        collectors = 1
    elif isinstance(sel, Selector):
        collectors = 1
    elif isinstance(sel, list):
        collectors = len(sel)
    if collectors == n:
        return True
    return False


class FmtType(Enum):
    REFERENCE    = auto()
    WILDCARD     = auto()
    FILENAME_REF = auto()
    FILENAME_GEN = auto()
    COMPLEX      = auto()

def get_fmttype(out: ToolOutput, tool: CommandTool) -> FmtType:
    # output uses WildcardSelector
    if isinstance(out.selector, WildcardSelector):
        return FmtType.WILDCARD

    # output uses InputSelector
    elif isinstance(out.selector, InputSelector):
        tinput = tool.inputs_map()[out.selector.input_to_select]
        
        # ToolInput is Filename type
        if isinstance(tinput.intype, Filename):
            entity_counts = trace_janis_entities(tinput.intype, tool=tool)
            entities = set(entity_counts.keys())
            filename_gen_whitelist = set(['Filename', 'str', 'NoneType'])
            filename_ref_whitelist = set(['InputSelector', 'Filename', 'str', 'NoneType'])
        
            # ToolInput does not refer to another ToolInput
            # This must be first as less specific
            if entities.issubset(filename_gen_whitelist):
                return FmtType.FILENAME_GEN
            
            # ToolInput refers to another ToolInput
            elif entities.issubset(filename_ref_whitelist):
                return FmtType.FILENAME_REF
            
            # ToolInput uses complex logic
            elif isinstance(tinput.intype, Filename):
                return FmtType.COMPLEX
        
        # ToolInput is not Filename type (direct reference)
        else:
            return FmtType.REFERENCE
    
    # anything else
    else:
        return FmtType.COMPLEX


### CMDTOOL OUTPUTS ###
class CmdtoolProcessOutputFactory:
    def __init__(self, out: ToolOutput, tool: CommandTool, sources: dict[str, Any]) -> None:
        self.out = out
        self.tool = tool
        self.sources = sources
        self.process_inputs = inputs.get_process_inputs(self.sources)
        self.param_inputs = inputs.get_param_inputs(self.sources)
        self.internal_inputs = inputs.get_internal_inputs(self.tool, self.sources)
        self.otype = get_otype(self.out)
        self.ftype = get_fmttype(self.out, self.tool)
        self.strategy_map = {
            OType.STDOUT: self.stdout_output,
            OType.NON_FILE: self.non_file_output,
            OType.FILE: self.file_output,
            OType.FILE_PAIR: self.file_pair_output,
            OType.FILE_ARRAY: self.file_array_output,
            OType.SECONDARIES: self.secondaries_output,
            OType.SECONDARIES_ARRAY: self.secondaries_array_output,
        }
        
        self.add_braces: bool = False
        self.add_quotes: bool = False


    # helper properties
    @property
    def basetype(self) -> DataType:
        return nfgen_utils.get_base_type(self.out.output_type)
    
    @property
    def dtype(self) -> DataType:
        return self.out.output_type
    
    @property
    def optional(self) -> bool:
        if self.dtype.optional:
            return True
        return False

    # helper methods
    def should_discard(self) -> bool:
        # TODO?
        return False
        
    def unwrap_collection_expression(self, expr: Any) -> str:
        # edge case - if referencing input (via InputSelector) and that
        # input is a Filename type, unwrap the Filename, not the InputSelector directly. 
        # results in cleaner format. 
        if isinstance(expr, InputSelector):
            inp = self.tool.inputs_map()[expr.input_to_select]
            if isinstance(inp.intype, Filename):
                expr = inp.intype

        if self.ftype == FmtType.REFERENCE:
            self.add_braces = False
            self.add_quotes = False
            expr = self.unwrap(expr)  
        elif self.ftype == FmtType.WILDCARD:
            self.add_braces = False
            self.add_quotes = False
            expr = self.unwrap(expr)  
            expr = f'"{expr}"'
        elif self.ftype == FmtType.FILENAME_GEN:
            self.add_braces = False
            self.add_quotes = False
            expr = self.unwrap(expr)
            expr = f'"{expr}"'
        elif self.ftype == FmtType.FILENAME_REF:
            self.add_braces = True
            self.add_quotes = False
            expr = self.unwrap(expr)
            expr = f'"{expr}"'
        elif self.ftype == FmtType.COMPLEX:
            self.add_braces = True
            self.add_quotes = False
            expr = self.unwrap(expr)
            expr = f'"{expr}"'
        else:
            raise NotImplementedError
        return expr

    def unwrap(self, expr: Any):
        return unwrap_expression(
            val=expr,
            tool=self.tool,
            sources=self.sources,
            process_inputs=self.process_inputs,
            param_inputs=self.param_inputs,
            internal_inputs=self.internal_inputs,
            in_shell_script=self.add_braces,
            quote_strings=self.add_quotes,
        )
    
    # public method
    def create(self) -> ProcessOutput:
        strategy = self.strategy_map[self.otype]
        process_output = strategy()
        return process_output

    # custom process output creation methods
    def stdout_output(self) -> StdoutProcessOutput:
        return StdoutProcessOutput(name=self.out.id(), is_optional=self.optional)
    
    def non_file_output(self) -> ValProcessOutput:
        expr = self.unwrap_collection_expression(self.out.selector)
        new_output = ValProcessOutput(
            name=self.out.id(), 
            is_optional=self.optional, 
            expression=expr
        )
        return new_output
    
    def file_output(self) -> PathProcessOutput:
        expr = self.unwrap_collection_expression(self.out.selector)
        new_output = PathProcessOutput(
            name=self.out.id(), 
            is_optional=self.optional, 
            expression=expr
        )
        return new_output
    
    def file_pair_output(self) -> TupleProcessOutput:
        assert(len(self.out.selector) == 2)
        qualifiers: list[str] = ['path', 'path']
        expressions: list[str] = []
        
        for item in self.out.selector:
            expr = self.unwrap_collection_expression(item)
            expressions.append(expr)
        
        new_output = TupleProcessOutput(
            name=self.out.id(), 
            is_optional=self.optional,
            qualifiers=qualifiers, 
            expressions=expressions
        )
        return new_output
    
    def file_array_output(self) -> PathProcessOutput:
        return self.file_output()
    
    def secondaries_output(self) -> TupleProcessOutput:
        """
        secondaries
        eg BamBai:
            selector=WildcardSelector("*.bam"),
            secondaries_present_as={".bai": ".bai"},
        """
        qualifiers: list[str] = []
        expressions: list[str] = []
        
        primary_expr = self.unwrap_collection_expression(self.out.selector)
        primary_expr = primary_expr.strip('"')
        exts = secondaries.get_extensions(self.dtype)
        for ext in exts:
            # primary file
            if self.out.secondaries_present_as is None or ext not in self.out.secondaries_present_as:
                qual = 'path'
                expr = f'"{primary_expr}"'
            # secondary file
            else:
                secondary_ext = self.out.secondaries_present_as[ext]
                secondary_expr: str = apply_secondary_file_format_to_filename(primary_expr, secondary_ext)
                qual = 'path'
                expr = f'"{secondary_expr}"'
            qualifiers.append(qual)
            expressions.append(expr)

        new_output = TupleProcessOutput(
            name=self.out.id(), 
            is_optional=self.optional,
            qualifiers=qualifiers, 
            expressions=expressions
        )
        return new_output
    
    def secondaries_array_output(self) -> SecondariesArrayProcessOutput:
        raise NotImplementedError
