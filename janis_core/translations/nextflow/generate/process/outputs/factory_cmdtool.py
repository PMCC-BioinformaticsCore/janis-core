

from typing import Any
from enum import Enum, auto

from janis_core.types import DataType
from janis_core.utils.secondary import apply_secondary_file_format_to_filename
from janis_core import (
    ToolOutput, 
    CommandTool,
    WildcardSelector,
    InputSelector,
    Selector,
    Filename
)
from janis_core import translation_utils as utils
from janis_core.translation_utils import DTypeType

from .... import task_inputs
from ....task_inputs import TaskInputType
from ....trace import trace_entity_counts
from ....unwrap import unwrap_expression
from ....variables import VariableManager

from ....model.process.outputs import (
    NFProcessOutput,
    NFPathProcessOutput,
    NFValProcessOutput,
    NFTupleProcessOutput,
    NFStdoutProcessOutput,
    NFSecondariesArrayProcessOutput
)

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
    FILENAME     = auto()  # filename TInput
    FILENAME_REF = auto()  # filename referencing another TInput: process input or param
    FILENAME_GEN = auto()  # filename referencing another TInput: internal input
    STATIC       = auto()  # static value for TInput
    COMPLEX      = auto()  # complex use of selectors / operators



### CMDTOOL OUTPUTS ###
class CmdtoolProcessOutputFactory:
    def __init__(
        self, 
        out: ToolOutput, 
        tool: CommandTool, 
        variable_manager: VariableManager,
    ) -> None:
        self.out = out
        self.tool = tool
        self.variable_manager = variable_manager
        self.dtt = utils.get_dtt(self.out.output_type)
        self.ftype = self.get_fmttype()
        self.strategy_map = {
            DTypeType.STDOUT: self.stdout_output,
            DTypeType.SECONDARY_ARRAY: self.secondaries_array_output,
            DTypeType.SECONDARY: self.secondaries_output,
            DTypeType.FILE_PAIR_ARRAY: self.file_pair_array_output,
            DTypeType.FILE_PAIR: self.file_pair_output,
            DTypeType.FILE_ARRAY: self.file_array_output,
            DTypeType.FILE: self.file_output,
            DTypeType.GENERIC_ARRAY: self.non_file_output,
            DTypeType.GENERIC: self.non_file_output,
        }
        
        self.add_braces: bool = False
        self.quote_strings: bool = False
    
    # public method
    def create(self) -> NFProcessOutput:
        strategy = self.strategy_map[self.dtt]
        process_output = strategy()
        return process_output
    
    # private below
    # helper properties
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
        if isinstance(self.out.selector, str):
            return FmtType.WILDCARD
        
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
                    if tinput.id() in task_inputs.task_inputs(self.tool.id()):
                        return FmtType.FILENAME
                    elif tinput.id() in task_inputs.param_inputs(self.tool.id()):
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
                # TInput with static value for process
                task_input = task_inputs.get(self.tool.id(), tinput)
                if task_input.ti_type == TaskInputType.STATIC:
                    return FmtType.STATIC
                # TInput which is referenced by variable
                else:
                    return FmtType.REFERENCE
        
        # anything else
        else:
            return FmtType.COMPLEX

    def unwrap_collection_expression(self, expr: Any) -> str:
        if self.ftype == FmtType.REFERENCE:
            self.add_braces = False
            self.quote_strings = False
            expr = self.unwrap(expr)  
        
        elif self.ftype in (FmtType.WILDCARD, FmtType.STATIC, FmtType.FILENAME_GEN):
            self.add_braces = False
            self.quote_strings = True
            expr = self.unwrap(expr)
            if not expr.startswith('"') and not expr.endswith('"'):
                expr = f'"{expr}"'
        
        elif self.ftype in (FmtType.FILENAME, FmtType.FILENAME_REF, FmtType.COMPLEX):
            self.add_braces = True
            self.quote_strings = True
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
            context='process_output',
            variable_manager=self.variable_manager,
            tool=self.tool,
            in_shell_script=self.add_braces,
            quote_strings=self.quote_strings,
        )
    
    # process output creation methods
    def stdout_output(self) -> NFStdoutProcessOutput:
        return NFStdoutProcessOutput(name=self.out.id(), janis_tag=self.out.id(), is_optional=self.optional)
    
    def non_file_output(self) -> NFValProcessOutput:
        expr = self.unwrap_collection_expression(self.out.selector)
        new_output = NFValProcessOutput(
            name=self.out.id(), 
            janis_tag=self.out.id(),
            is_optional=self.optional, 
            expression=expr
        )
        return new_output
    
    def file_output(self) -> NFPathProcessOutput:
        expr = self.unwrap_collection_expression(self.out.selector)
        new_output = NFPathProcessOutput(
            name=self.out.id(), 
            janis_tag=self.out.id(),
            is_optional=self.optional, 
            expression=expr
        )
        return new_output
    
    def file_pair_output(self) -> NFPathProcessOutput | NFTupleProcessOutput:
        if self.ftype == FmtType.WILDCARD:
            return self.file_output()
        elif has_n_collectors(self.out, n=2):
            expr1 = self.unwrap_collection_expression(self.out.selector[0])
            expr2 = self.unwrap_collection_expression(self.out.selector[1])
            new_output = NFTupleProcessOutput(
                name=self.out.id(), 
                janis_tag=self.out.id(),
                is_optional=self.optional,
                qualifiers=['path', 'path'], 
                expressions=[expr1, expr2]
            )
            return new_output
        else:
            raise NotImplementedError
    
    def file_pair_array_output(self) -> NFPathProcessOutput:
        if self.ftype == FmtType.WILDCARD:
            # TODO raise warning
            return self.file_output()
        else:
            raise NotImplementedError
        
    def file_array_output(self) -> NFPathProcessOutput:
        return self.file_output()
    
    def secondaries_output(self) -> NFTupleProcessOutput:
        """
        eg BamBai:
            selector=WildcardSelector("*.bam"),
            secondaries_present_as={".bai": ".bai"},
        """
        if self.ftype == FmtType.REFERENCE:
            return self.secondaries_output_reference()
        else:
            return self.secondaries_output_complex()

    def secondaries_output_reference(self) -> NFTupleProcessOutput:
        qualifiers: list[str] = []
        expressions: list[str] = []
        
        # primary file
        primary_reference = self.unwrap_collection_expression(self.out.selector)
        qualifiers.append('path')
        expressions.append(primary_reference)

        exts = utils.get_extensions(self.dtype)               
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

        new_output = NFTupleProcessOutput(
            name=self.out.id(), 
            janis_tag=self.out.id(),
            is_optional=self.optional,
            qualifiers=qualifiers, 
            expressions=expressions
        )
        return new_output
    
    def secondaries_output_complex(self) -> NFTupleProcessOutput:
        qualifiers: list[str] = []
        expressions: list[str] = []

        exts = utils.get_extensions(self.dtype)
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

        new_output = NFTupleProcessOutput(
            name=self.out.id(), 
            janis_tag=self.out.id(),
            is_optional=self.optional,
            qualifiers=qualifiers, 
            expressions=expressions
        )
        return new_output

    def secondaries_array_output(self) -> NFSecondariesArrayProcessOutput:
        raise NotImplementedError('process outputs with format [[file1, file2]] (arrays of secondary files) not supported in nextflow translation')
