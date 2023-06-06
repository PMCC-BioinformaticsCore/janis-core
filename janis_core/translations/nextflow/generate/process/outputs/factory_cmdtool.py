

import regex as re
from typing import Any, Optional
from enum import Enum, auto

from janis_core.types import DataType
from janis_core import (
    ToolOutput, 
    CommandTool,
    WildcardSelector,
    InputSelector,
    Selector,
    Filename,
    Stdout
)
from janis_core import translation_utils as utils
from janis_core.translation_utils import DTypeType

from janis_core.translations.common import trace
from .... import task_inputs
from ....task_inputs import TaskInputType
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
        self.strategy_map = {
            DTypeType.STDOUT: self.stdout_output,
            DTypeType.SECONDARY_ARRAY: self.secondaries_array_output,
            DTypeType.SECONDARY: self.secondaries_output_dispatcher,
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
    def get_fmttype(self, selector: Any) -> FmtType:
        """returns a FmtType based on the specific ToolOutput we have received"""
        # output uses WildcardSelector
        if isinstance(selector, str):
            return FmtType.WILDCARD
        
        if isinstance(selector, WildcardSelector):
            return FmtType.WILDCARD

        # output uses InputSelector
        elif isinstance(selector, InputSelector):
            tinput = self.tool.inputs_map()[selector.input_to_select]
            
            # ToolInput is Filename type
            if isinstance(tinput.intype, Filename):
                entity_counts = trace.trace_entity_counts(tinput.intype, tool=self.tool)
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

    def unwrap_collection_expression(self, expr: Any, ftype: FmtType) -> str:
        if ftype == FmtType.REFERENCE:
            self.add_braces = False
            self.quote_strings = False
            expr = self.unwrap(expr)  
        
        elif ftype in (FmtType.WILDCARD, FmtType.STATIC, FmtType.FILENAME_GEN):
            self.add_braces = False
            self.quote_strings = True
            expr = self.unwrap(expr)
            if not expr.startswith('"') and not expr.endswith('"'):
                expr = f'"{expr}"'
        
        elif ftype in (FmtType.FILENAME, FmtType.FILENAME_REF, FmtType.COMPLEX):
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
    def stdout_output(self) -> NFProcessOutput:
        dtt = utils.get_dtt(self.out.output_type.subtype)
        if dtt in [
            DTypeType.GENERIC_ARRAY,
            DTypeType.GENERIC,
        ]:
            return NFStdoutProcessOutput(
                name=self.out.id(), 
                janis_tag=self.out.id(), 
                is_optional=self.optional
            )
        else:
            if hasattr(self.out.output_type.subtype, 'extension') and self.out.output_type.subtype.extension is not None:
                suffix = self.out.output_type.subtype.extension
            else:
                suffix = ''
            return NFPathProcessOutput(
                name=self.out.id(), 
                janis_tag=self.out.id(), 
                is_optional=self.optional,
                expression=f'"{self.out.id()}{suffix}"'
            )
    
    def non_file_output(self) -> NFValProcessOutput:
        ftype = self.get_fmttype(self.out.selector)
        expr = self.unwrap_collection_expression(self.out.selector, ftype)
        new_output = NFValProcessOutput(
            name=self.out.id(), 
            janis_tag=self.out.id(),
            is_optional=self.optional, 
            expression=expr
        )
        return new_output
    
    def file_output(self) -> NFPathProcessOutput:
        ftype = self.get_fmttype(self.out.selector)
        expr = self.unwrap_collection_expression(self.out.selector, ftype)
        new_output = NFPathProcessOutput(
            name=self.out.id(), 
            janis_tag=self.out.id(),
            is_optional=self.optional, 
            expression=expr
        )
        return new_output
    
    def file_pair_output(self) -> NFPathProcessOutput | NFTupleProcessOutput:
        ftype = self.get_fmttype(self.out.selector)
        if ftype == FmtType.WILDCARD:
            return self.file_output()
        elif has_n_collectors(self.out, n=2):
            ftype1 = self.get_fmttype(self.out.selector[0])
            ftype2 = self.get_fmttype(self.out.selector[1])
            expr1 = self.unwrap_collection_expression(self.out.selector[0], ftype1)
            expr2 = self.unwrap_collection_expression(self.out.selector[1], ftype2)
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
        ftype = self.get_fmttype(self.out.selector)
        if ftype == FmtType.WILDCARD:
            # TODO raise warning
            return self.file_output()
        else:
            raise NotImplementedError
        
    def file_array_output(self) -> NFPathProcessOutput:
        return self.file_output()
    
    def secondaries_output_dispatcher(self) -> NFTupleProcessOutput:
        """
        eg BamBai:
            selector=WildcardSelector("*.bam"),
            secondaries_present_as={".bai": ".bai"},
        """
        exts = utils.get_extensions(self.dtype, ignore_duplicates=True)
        
        if has_n_collectors(self.out, n=1):
            return self.secondaries_output_single_collector()
        
        elif has_n_collectors(self.out, n=len(exts)):
            return self.secondaries_output_n_collectors()
        
        else:
            raise NotImplementedError

    def secondaries_output_single_collector(self) -> NFTupleProcessOutput:
        exts = utils.get_extensions(self.dtype, ignore_duplicates=True)
        
        qualifiers: list[str] = ['path'] * len(exts)
        expressions: list[str] = []

        base_expr = self.resolve_secondary_collector(self.out.selector)
        expressions.append(base_expr)
        
        for ext in exts[1:]:
            expr = self.resolve_secondary_collector(self.out.selector, ext)
            expressions.append(expr)

        new_output = NFTupleProcessOutput(
            name=self.out.id(), 
            janis_tag=self.out.id(),
            is_optional=self.optional,
            qualifiers=qualifiers, 
            expressions=expressions
        )
        return new_output
    
    def secondaries_output_n_collectors(self) -> NFTupleProcessOutput:
        qualifiers: list[str] = ['path'] * len(self.out.selector)
        expressions: list[str] = []

        for sel in self.out.selector:
            expr = self.resolve_secondary_collector(sel)
            expressions.append(expr)
        
        new_output = NFTupleProcessOutput(
            name=self.out.id(), 
            janis_tag=self.out.id(),
            is_optional=self.optional,
            qualifiers=qualifiers, 
            expressions=expressions
        )
        return new_output
    
    def resolve_secondary_collector(self, entity: Any, ext: Optional[str]=None) -> str:
        ftype = self.get_fmttype(entity)
        expr = self.unwrap_collection_expression(entity, ftype)
        if ext is not None:
            try:
                expr = self.replace_extension_in_expression_unsafe(expr, ext)
            except Exception as e:
                print(expr)
                print(e)
                expr = self.replace_extension_in_expression_safe(expr, ext)

        return expr
    
    def replace_extension_in_expression_unsafe(self, expr: str, ext: str) -> str:
        if not expr.startswith('"') and not expr.endswith('"'):
            expr = f'"{expr}"'
        matcher = r"\$\{(.+)\}"
        matches = list(re.finditer(matcher, expr))
        if len(matches) == 1:
            match = matches[0]
            if ext.startswith('^') and not match.groups()[0].endswith('.simpleName'):
                ext = ext.lstrip('^')
                bracket_contents_start = match.regs[1][0]
                bracket_contents = match.groups()[0]
                expr = expr[:bracket_contents_start] + bracket_contents + '.simpleName' + f'}}{ext}' + expr[match.end():]
                print()
            elif ext.startswith('^'):
                ext = ext.lstrip('^')
                expr = expr[:match.end()] + ext + expr[match.end():]
                print()
            else:
                expr = expr[:match.end()] + ext + expr[match.end():]
            return expr
        else:
            raise ValueError()

    def replace_extension_in_expression_safe(self, expr: str, ext: str) -> str:
        ext = ext.lstrip('^')
        return f'"*{ext}"'

    def secondaries_array_output(self) -> NFSecondariesArrayProcessOutput:
        raise NotImplementedError('process outputs with format [[file1, file2]] (arrays of secondary files) not supported in nextflow translation')

