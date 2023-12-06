

import regex as re
from typing import Any, Optional
from enum import Enum, auto

from janis_core.types import DataType
from janis_core import (
    ToolOutput, 
    CommandToolBuilder,
    WildcardSelector,
    InputSelector,
    Selector,
    Filename,
    Stdout,
    StringFormatter,
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
    STRINGFORMATTER = auto()  # string formatter
    REFERENCE       = auto()  # reference to process input or param
    WILDCARD        = auto()  # regex based collection
    FILENAME        = auto()  # filename TInput
    STATIC          = auto()  # static value for TInput
    COMPLEX         = auto()  # complex use of selectors / operators


### CMDTOOL OUTPUTS ###
class CmdtoolProcessOutputFactory:
    def __init__(
        self, 
        out: ToolOutput, 
        tool: CommandToolBuilder, 
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
            DTypeType.FILENAME: self.non_file_output,
            DTypeType.GENERIC_ARRAY: self.non_file_output,
            DTypeType.GENERIC: self.non_file_output,
        }
    
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
        
        elif isinstance(selector, WildcardSelector):
            return FmtType.WILDCARD
        
        elif isinstance(selector, StringFormatter):
            return FmtType.STRINGFORMATTER

        # output uses InputSelector
        elif isinstance(selector, InputSelector):
            if selector.input_to_select not in self.tool.inputs_map():
                print()
            tinput = self.tool.inputs_map()[selector.input_to_select]
            
            # ToolInput is Filename type
            if isinstance(tinput.intype, Filename):
                return FmtType.FILENAME 
                
            # ToolInput is not Filename type (direct reference)
            else:
                # TInput with static value for process
                task_input = task_inputs.get(self.tool.id(), tinput)
                if task_input.ti_type == TaskInputType.STATIC:
                    return FmtType.STATIC
                # TInput which is directly referenced by variable
                else:
                    return FmtType.REFERENCE
        
        # anything else
        else:
            return FmtType.COMPLEX

    def unwrap_collection_expression(self, expr: Any, ftype: FmtType) -> str:
        if ftype == FmtType.REFERENCE:
            expr = self.unwrap(expr) 
        
        elif ftype == FmtType.STRINGFORMATTER:
            expr = self.unwrap(expr, apply_braces=False) 
            expr = f'"{expr}"'
        
        elif ftype in (FmtType.WILDCARD, FmtType.STATIC):
            expr = self.unwrap(expr)
            if not expr.startswith('"') and not expr.endswith('"'):
                expr = f'"{expr}"'
        
        # things which need unwrapping
        elif ftype in (FmtType.FILENAME, FmtType.COMPLEX):
            expr = self.unwrap(expr, apply_braces=True)
            if not expr.startswith('"') and not expr.endswith('"'):
                expr = f'"{expr}"'
        
        else:
            # some future FmtType
            raise NotImplementedError
        return expr

    def unwrap(self, expr: Any, apply_braces: bool=False):
        return unwrap_expression(
            val=expr,
            context='process_output',
            variable_manager=self.variable_manager,
            tool=self.tool,
            apply_braces=apply_braces
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
        eg BamBai replaced:
            selector=WildcardSelector("*.bam"),
            secondaries_present_as={".bai": "^.bai"},
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
        if self.out.secondaries_present_as:
            exts = [exts[0]] + [self.out.secondaries_present_as[x] for x in exts[1:]]
        
        qualifiers: list[str] = ['path'] * len(exts)
        expressions: list[str] = []

        base_expr = self.resolve_secondary_collector(self.out.selector)
        expressions.append(base_expr)
        
        for ext in exts[1:]:
            expr = self.resolve_secondary_collector(self.out.selector, primary_ext=exts[0], secondary_ext=ext)
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
    
    def resolve_secondary_collector(self, entity: Any, primary_ext: Optional[str]=None, secondary_ext: Optional[str]=None) -> str:
        ftype = self.get_fmttype(entity)
        expr = self.unwrap_collection_expression(entity, ftype)
        if secondary_ext is not None:
            assert(primary_ext is not None)
            try:
                expr = self.replace_extension_in_expression_unsafe(expr, primary_ext, secondary_ext)
            except Exception as e:
                expr = self.replace_extension_in_expression_safe(expr, primary_ext, secondary_ext)
        return expr
    
    def replace_extension_in_expression_unsafe(self, expr: str, primary_ext: str, secondary_ext: str) -> str:
        primary_ext = primary_ext.replace('.', '\.')
        primary_ext_matcher = fr"{primary_ext}(?=['\"\s])"
        matches = list(re.finditer(primary_ext_matcher, expr))
        if len(matches) == 1:
            match = matches[0]
            if secondary_ext.startswith('^'):
                expr = expr[:match.start()] + secondary_ext.lstrip('^') + expr[match.end():]
            else:
                expr = expr[:match.end()] + secondary_ext + expr[match.end():]
        else:
            raise ValueError()
        return expr

    def replace_extension_in_expression_safe(self, expr: str, primary_ext: str, secondary_ext: str) -> str:
        ext = secondary_ext.lstrip('^')
        return f'"*{ext}"'

    def secondaries_array_output(self) -> NFSecondariesArrayProcessOutput:
        raise NotImplementedError('process outputs with format [[file1, file2]] (arrays of secondary files) not supported in nextflow translation')

