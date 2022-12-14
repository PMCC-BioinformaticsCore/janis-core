
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any
# from functools import cached_property

from janis_core.types import File, Directory, Array, DataType, Stdout
from janis_core import (
    ToolOutput, 
    TOutput, 
    CommandTool, 
    PythonTool,
)

from janis_core.utils.secondary import apply_secondary_file_format_to_filename

from ..unwrap import unwrap_expression
from .. import nfgen_utils
from .. import secondaries
from .. import settings
from . import inputs


@dataclass
class ProcessOutput(ABC):
    name: str
    is_optional: bool

    @abstractmethod
    def get_string(self) -> str:
        ...
    
    @property
    def emit(self) -> str:
        return f', emit: {self.name}'
    
    @property
    def optional(self) -> str:
        if self.is_optional:
            return ', optional: true'
        return ''


@dataclass
class StdoutProcessOutput(ProcessOutput):

    def get_string(self) -> str:
        return f'stdout{self.optional}{self.emit}'


@dataclass
class ValProcessOutput(ProcessOutput):
    expression: str

    def get_string(self) -> str:
        return f'val {self.expression}{self.optional}{self.emit}'


@dataclass
class PathProcessOutput(ProcessOutput):
    expression: str

    def get_string(self) -> str:
        return f'path {self.expression}{self.optional}{self.emit}'


@dataclass
class TupleProcessOutput(ProcessOutput):
    qualifiers: list[str]
    expressions: list[str]

    @property
    def fields(self) -> str:
        out: str = ''
        for qual, expr in zip(self.qualifiers, self.expressions):
            out += f'{qual}({expr}), '
        out = out.rstrip(', ') # strip off the last comma & space
        return out
    
    def get_string(self) -> str:
        return f'tuple {self.fields}{self.optional}{self.emit}'

 

# def filter_null(process_outs: list[Optional[ProcessOutput]]) -> list[ProcessOutput]:
#     return [x for x in process_outs if x is not None]



def create_outputs(out: ToolOutput | TOutput, tool: CommandTool | PythonTool, sources: dict[str, Any]) -> list[ProcessOutput]:
    if isinstance(out, ToolOutput) and isinstance(tool, CommandTool):
        factory = CmdtoolProcessOutputFactory(out, tool, sources)
    if isinstance(out, TOutput) and isinstance(tool, PythonTool):
        factory = PythonToolProcessOutputFactory(out, tool, sources)
    return factory.create()


### CMDTOOL OUTPUTS ###
class CmdtoolProcessOutputFactory:
    def __init__(self, out: ToolOutput, tool: CommandTool, sources: dict[str, Any]) -> None:
        self.out = out
        self.tool = tool
        self.sources = sources
        self.process_inputs = inputs.get_process_inputs(self.sources)
        self.param_inputs = inputs.get_param_inputs(self.sources)
        self.internal_inputs = inputs.get_internal_inputs(self.tool, self.sources)

    @property
    def basetype(self) -> DataType:
        return nfgen_utils.get_base_type(self.out.output_type)
    
    @property
    def dtype(self) -> DataType:
        return self.out.output_type

    def unwrap_collection_expression(self) -> Any:
        return unwrap_expression(
            val=self.out.selector, 
            tool=self.tool, 
            in_shell_script=True,
            sources=self.sources,
            process_inputs=self.process_inputs,
            param_inputs=self.param_inputs,
            internal_inputs=self.internal_inputs,
        )

    def create(self) -> list[ProcessOutput]:
        if self.should_discard():
            return []
        if isinstance(self.dtype, Array):
            return self.create_outputs_array()
        else:
            return self.create_outputs_single()

    def should_discard(self) -> bool:
        # TODO?
        return False
    
    def create_outputs_single(self) -> list[ProcessOutput]:
        process_outs: list[ProcessOutput] = []
        # stdout
        if isinstance(self.basetype, Stdout):
            process_outs = [self.create_stdout_output()]

        # file secondaries
        elif secondaries.is_secondary_type(self.dtype):
            process_outs = [self.create_tuple_output_secondaries()]
        
        # file
        elif isinstance(self.basetype, (File, Directory)):
            process_outs = [self.create_path_output()]
        
        # nonfile
        else:
            process_outs = [self.create_val_output()]
        
        return process_outs

    def create_outputs_array(self) -> list[ProcessOutput]:
        process_outs: list[ProcessOutput] = []

        # secondaries array
        if secondaries.is_array_secondary_type(self.dtype):
            # a path output per file type
            assert(isinstance(self.basetype, File))
            exts = secondaries.get_names(self.basetype)
            for ext in exts:
                process_outs.append(self.create_path_output_secondaries(ext))

        # file array
        elif isinstance(self.basetype, (File, Directory)):
            process_outs = [self.create_path_output()]

        # nonfile array
        else:
            process_outs = [self.create_val_output()]

        return process_outs

    def create_stdout_output(self) -> StdoutProcessOutput:
        # stdout
        optional = True if self.dtype.optional else False  
        new_output = StdoutProcessOutput(name=self.out.id(), is_optional=optional)
        return new_output

    def create_path_output(self) -> PathProcessOutput:
        # create output for file types
        optional = True if self.dtype.optional else False
        expression = self.unwrap_collection_expression()
        expression = f'"{expression}"'
        new_output = PathProcessOutput(
            name=self.out.id(), 
            is_optional=optional, 
            expression=expression
        )
        return new_output

    def create_val_output(self) -> ValProcessOutput:
        # create output for nonfile types
        optional = True if self.dtype.optional else False  
        expression = self.unwrap_collection_expression()
        expression = f'"{expression}"'
        new_output = ValProcessOutput(
            name=self.out.id(), 
            is_optional=optional, 
            expression=expression
        )
        return new_output

    def create_tuple_output_secondaries(self) -> TupleProcessOutput:
        """
        secondaries
        eg BamBai:
            selector=WildcardSelector("*.bam"),
            secondaries_present_as={".bai": ".bai"},
        """
        assert(isinstance(self.dtype, File))
        optional = True if self.dtype.optional else False  
        qualifiers: list[str] = []
        expressions: list[str] = []
        
        primary_expr = self.unwrap_collection_expression()
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
            is_optional=optional,
            qualifiers=qualifiers, 
            expressions=expressions
        )
        return new_output

    def create_path_output_secondaries(self, ext: str) -> PathProcessOutput:
        # array of secondaries
        raise NotImplementedError




"""
HMMMMM

    # def should_discard(self) -> bool:
    #     if self.out.selector and isinstance(self.out.selector, InputSelector):
    #         if 
    #     return False


        # if self.out.selector is not None and isinstance(self.out.selector, WildcardSelector):
        #     expression = self.get_wildcard_selector_str()
        
        # elif self.out.selector is not None and isinstance(self.out.selector, InputSelector):
        #     expression = self.get_input_selector_str()
        
        # else:
        #     raise NotImplementedError

    # def get_wildcard_selector_str(self) -> str:
    #     selector: WildcardSelector = self.out.selector # type: ignore
    #     return f'"{selector.wildcard}"'
    
    # def get_input_selector_str(self) -> str:
    #     expression = ''

    #     selector: InputSelector = self.out.selector # type: ignore
    #     inname = selector.input_to_select
    #     if inname in self.process_inputs:
    #         expression = inname
            
    #     elif inname in self.param_inputs:
    #         # get the param?
    #         src = self.sources[inname]
    #         sel = src.source_map[0].source
    #         raise NotImplementedError
    #         param = params.get(sel.input_node.uuid)
    #         return f'params.{param.name}'
        
    #     elif inname in self.internal_inputs:
    #         pass
        
    #     else:
    #         raise NotImplementedError

    #     # remove_file_extension
    #     if selector.remove_file_extension:
    #         expression = f"{expression}.simpleName"

"""


### PYTHONTOOL OUTPUTS ###
class PythonToolProcessOutputFactory:
    def __init__(self, out: TOutput, tool: PythonTool, sources: dict[str, Any]) -> None:
        self.out = out
        self.tool = tool
        self.sources = sources
        self.process_inputs = inputs.get_process_inputs(self.sources)
        self.param_inputs = inputs.get_param_inputs(self.sources)
        self.internal_inputs = inputs.get_internal_inputs(self.tool, self.sources)

    @property
    def basetype(self) -> DataType:
        return nfgen_utils.get_base_type(self.out.outtype)
    
    @property
    def dtype(self) -> DataType:
        return self.out.outtype

    def create(self) -> list[ProcessOutput]:
        if isinstance(self.dtype, Array):
            return self.create_outputs_array()
        else:
            return self.create_outputs_single()
    
    def create_outputs_single(self) -> list[ProcessOutput]:
        # file secondaries
        if secondaries.is_secondary_type(self.dtype):
            raise NotImplementedError
            # outputs = [create_tuple_output_secondaries(out, tool)]
            # return outputs
        
        # file
        if isinstance(self.basetype, (File, Directory)):
            return [self.create_path_output()]
        
        # nonfile
        return [self.create_val_output()]

    def create_outputs_array(self) -> list[ProcessOutput]:
        # secondaries array
        if secondaries.is_array_secondary_type(self.dtype):
            # a path output per file type
            outputs: list[ProcessOutput] = []
            assert(isinstance(self.basetype, File))
            exts = secondaries.get_names(self.basetype)
            for ext in exts:
                outputs.append(self.create_path_output_secondaries(ext))
            return outputs

        # file array
        if isinstance(self.basetype, (File, Directory)):
            return [self.create_path_output()]

        # nonfile array
        return [self.create_val_output_array()]

    def create_path_output(self) -> PathProcessOutput:
        # file
        optional = True if self.dtype.optional else False  
        filename = f'"{settings.PYTHON_CODE_OUTPUT_FILENAME_PREFIX}{self.out.id()}"'
        new_output = PathProcessOutput(
            name=self.out.id(), 
            is_optional=optional, 
            expression=filename
        )
        return new_output

    def create_val_output(self) -> ValProcessOutput:
        # nonfile
        optional = True if self.dtype.optional else False  
        filename = f'{settings.PYTHON_CODE_OUTPUT_FILENAME_PREFIX}{self.out.id()}'
        expression = f'\"${{file(\"${{task.workDir}}/{filename}\").text}}"'
        new_output = ValProcessOutput(
            name=self.out.id(), 
            is_optional=optional, 
            expression=expression
        )
        return new_output

    def create_val_output_array(self) -> ValProcessOutput:
        # nonfile
        optional = True if self.dtype.optional else False  
        filename = f'{settings.PYTHON_CODE_OUTPUT_FILENAME_PREFIX}{self.out.id()}'
        expression = f'\"${{file(\"${{task.workDir}}/{filename}\").text.split(\',\')}}"'
        new_output = ValProcessOutput(
            name=self.out.id(), 
            is_optional=optional, 
            expression=expression
        )
        return new_output

    def create_path_output_secondaries(self, ext: str) -> PathProcessOutput:
        # array of secondaries
        raise NotImplementedError



