

from typing import Any
from enum import Enum, auto

from janis_core.types import File, DataType, Stdout, Directory
from janis_core.utils.secondary import apply_secondary_file_format_to_filename
from janis_core import (
    ToolOutput, 
    CommandTool,
    WildcardSelector,
    InputSelector,
    Selector,
    Filename
)

from ...plumbing import trace_entity_counts
from ...unwrap import unwrap_expression
from ... import nfgen_utils
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
    FILEPAIR            = auto()
    FILE_ARRAY          = auto()
    SECONDARIES         = auto()
    SECONDARIES_ARRAY   = auto()


def get_otype(out: ToolOutput) -> OType:
    if is_stdout_type(out):
        return OType.STDOUT

    elif is_filepair_type(out):
        return OType.FILEPAIR
    
    elif is_secondary_type(out) and is_array_type(out):
        return OType.SECONDARIES_ARRAY
    
    elif is_secondary_type(out):
        return OType.SECONDARIES
    
    elif is_file_type(out) and is_array_type(out) and has_n_collectors(out, n=1):
        return OType.FILE_ARRAY
    
    elif is_file_type(out):
        return OType.FILE
    
    elif is_non_file_type(out) and has_n_collectors(out, n=1):
        return OType.NON_FILE

    else:
        # some future OType
        raise NotImplementedError


def is_stdout_type(out: ToolOutput) -> bool:
    if isinstance(out.output_type, Stdout):
        return True
    return False

def is_file_type(out: ToolOutput) -> bool:
    basetype = nfgen_utils.get_base_type(out.output_type)
    basetype = nfgen_utils.ensure_single_type(basetype)
    if isinstance(basetype, (File, Directory)):
        return True
    return False

def is_filepair_type(out: ToolOutput) -> bool:
    basetype = nfgen_utils.get_base_type(out.output_type)
    basetype = nfgen_utils.ensure_single_type(basetype)
    if basetype.name() in ['FastqPair', 'FastqGzPair']:
        return True
    return False

def is_secondary_type(out: ToolOutput) -> bool:
    return nfgen_utils.is_secondary_type(out.output_type)

def is_array_type(out: ToolOutput) -> bool:
    if out.output_type.is_array():
        return True
    return False

def is_non_file_type(out: ToolOutput) -> bool:
    basetype = nfgen_utils.get_base_type(out.output_type)
    basetype = nfgen_utils.ensure_single_type(basetype)
    if not isinstance(basetype, (File, Directory)):
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
    REFERENCE    = auto()  # reference to process input or param
    WILDCARD     = auto()  # regex based collection
    FILENAME     = auto()  # filename ToolInput
    FILENAME_REF = auto()  # filename referencing another ToolInput: process input or param
    FILENAME_GEN = auto()  # filename referencing another ToolInput: internal input
    COMPLEX      = auto()  # complex use of selectors / operators



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
        self.ftype = self.get_fmttype()
        self.strategy_map = {
            OType.STDOUT: self.stdout_output,
            OType.NON_FILE: self.non_file_output,
            OType.FILE: self.file_output,
            OType.FILEPAIR: self.file_pair_output,
            OType.FILE_ARRAY: self.file_array_output,
            OType.SECONDARIES: self.secondaries_output,
            OType.SECONDARIES_ARRAY: self.secondaries_array_output,
        }
        
        self.add_braces: bool = False
        self.add_quotes: bool = False
    
    # public method
    def create(self) -> ProcessOutput:
        strategy = self.strategy_map[self.otype]
        process_output = strategy()
        return process_output
    
    # private below
    # helper properties
    @property
    def basetype(self) -> DataType:
        basetype = nfgen_utils.get_base_type(self.out.output_type)
        basetype = nfgen_utils.ensure_single_type(basetype)
        return basetype
    
    @property
    def dtype(self) -> DataType:
        return self.out.output_type
    
    @property
    def optional(self) -> bool:
        if self.dtype.optional:
            return True
        return False

    # helper methods
    def get_fmttype(self) -> FmtType:
        """returns a FmtType based on the specific ToolOutput we have received"""
        # output uses WildcardSelector
        if isinstance(self.out.selector, WildcardSelector):
            return FmtType.WILDCARD

        # output uses InputSelector
        elif isinstance(self.out.selector, InputSelector):
            tinput = self.tool.inputs_map()[self.out.selector.input_to_select]
            
            # ToolInput is Filename type
            if isinstance(tinput.intype, Filename):
                entity_counts = trace_entity_counts(tinput.intype, tool=self.tool)
                entities = set(entity_counts.keys())
                filename_gen_whitelist = set(['Filename', 'str', 'NoneType'])
                filename_ref_whitelist = set(['InputSelector', 'Filename', 'str', 'NoneType'])
            
                # ToolInput does not refer to another ToolInput
                # This must be first as less specific
                if entities.issubset(filename_gen_whitelist):
                    if tinput.id() in self.process_inputs or tinput.id() in self.param_inputs:
                        return FmtType.FILENAME
                    else:
                        return FmtType.FILENAME_GEN
                
                # ToolInput refers to another ToolInput
                elif entities.issubset(filename_ref_whitelist):
                    return FmtType.FILENAME_REF
                
                # ToolInput uses complex logic
                else:
                    return FmtType.COMPLEX
            
            # ToolInput is not Filename type (direct reference)
            else:
                return FmtType.REFERENCE
        
        # anything else
        else:
            return FmtType.COMPLEX

    def unwrap_collection_expression(self, expr: Any) -> str:
        if self.ftype == FmtType.REFERENCE:
            self.add_braces = False
            self.add_quotes = False
            expr = self.unwrap(expr)  
        elif self.ftype in (FmtType.WILDCARD, FmtType.FILENAME_GEN):
            self.add_braces = False
            self.add_quotes = False
            expr = self.unwrap(expr)
            if not expr.startswith('"') and not expr.endswith('"'):
                expr = f'"{expr}"'
        elif self.ftype in (FmtType.FILENAME, FmtType.FILENAME_REF, FmtType.COMPLEX):
            self.add_braces = True
            self.add_quotes = False
            expr = self.unwrap(expr)
            if not expr.startswith('"') and not expr.endswith('"'):
                expr = f'"{expr}"'
        else:
            # some future FmtType
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
    
    # process output creation methods
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
        eg BamBai:
            selector=WildcardSelector("*.bam"),
            secondaries_present_as={".bai": ".bai"},
        """
        if self.ftype == FmtType.REFERENCE:
            return self.secondaries_output_reference()
        else:
            return self.secondaries_output_complex()


    def secondaries_output_reference(self) -> TupleProcessOutput:
        qualifiers: list[str] = []
        expressions: list[str] = []
        
        # primary file
        primary_reference = self.unwrap_collection_expression(self.out.selector)
        qualifiers.append('path')
        expressions.append(primary_reference)

        exts = nfgen_utils.get_extensions(self.dtype)               
        for ext in exts[1:]:
            # extensions which replace part of filename
            if ext.startswith('^'):
                secondary_ext = ext.lstrip('^')
                qualifiers.append('path')
                expressions.append(f'"*{secondary_ext}"')
            else:
                secondary_ext = ext
                qualifiers.append('path')
                expressions.append(f'"${{{primary_reference}}}{secondary_ext}"')

        new_output = TupleProcessOutput(
            name=self.out.id(), 
            is_optional=self.optional,
            qualifiers=qualifiers, 
            expressions=expressions
        )
        return new_output
    
    def secondaries_output_complex(self) -> TupleProcessOutput:
        qualifiers: list[str] = []
        expressions: list[str] = []

        exts = nfgen_utils.get_extensions(self.dtype)
        primary_expr = self.unwrap_collection_expression(self.out.selector)
        primary_expr = primary_expr.strip('"')
        primary_ext = exts[0]
        primary_ext_start = primary_expr.rfind(primary_ext)
        primary_ext_stop = primary_ext_start + len(primary_ext)

        # primary file
        qualifiers.append('path')
        expressions.append(f'"{primary_expr}"')

        # secondary files
        for ext in exts[1:]:
            # no primary ext found in collection expression (edge case)
            if primary_ext_start == -1:
                # TODO add a user warning message - dodgy output collection in this step
                secondary_ext = ext.lstrip('^')
                qualifiers.append('path')
                expressions.append(f'"*{secondary_ext}"')

            # find the last occurance of the primary file format, replace this with secondary file format
            # eg "${alignment}.bam" -> "${alignment}.bai"
            # probably has bugs.
            else:
                if self.out.secondaries_present_as and ext in self.out.secondaries_present_as:
                    secondary_ext_pattern = self.out.secondaries_present_as[ext]
                else:
                    secondary_ext_pattern = ext
                secondary_ext: str = apply_secondary_file_format_to_filename(primary_ext, secondary_ext_pattern)
                secondary_expr = primary_expr[:primary_ext_start] + secondary_ext + primary_expr[primary_ext_stop:]
                qualifiers.append('path')
                expressions.append(f'"{secondary_expr}"')

        new_output = TupleProcessOutput(
            name=self.out.id(), 
            is_optional=self.optional,
            qualifiers=qualifiers, 
            expressions=expressions
        )
        return new_output

    def secondaries_array_output(self) -> SecondariesArrayProcessOutput:
        raise NotImplementedError('process outputs with format [[file1, file2]] (arrays of secondary files) not supported in nextflow translation')
