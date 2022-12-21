

from typing import Any
from enum import Enum, auto

from janis_core.types import File
from janis_core import (
    TOutput, 
    PythonTool,
)

from ... import settings
from ... import nfgen_utils

from .model import (
    ProcessOutput,
    PathProcessOutput,
    ValProcessOutput,
)




class OType(Enum):
    NON_FILE        = auto()
    NON_FILE_ARRAY  = auto()
    FILE            = auto()


def get_otype(out: TOutput) -> OType:
    if is_file_type(out):
        return OType.FILE
    elif is_non_file_type(out) and is_array_type(out):
        return OType.NON_FILE_ARRAY
    elif is_non_file_type(out):
        return OType.NON_FILE
    else:
        raise NotImplementedError

def is_file_type(out: TOutput) -> bool:
    basetype = nfgen_utils.get_base_type(out.outtype)
    if isinstance(basetype, File):
        return True
    return False

def is_array_type(out: TOutput) -> bool:
    if out.outtype.is_array():
        return True
    return False

def is_non_file_type(out: TOutput) -> bool:
    basetype = nfgen_utils.get_base_type(out.outtype)
    if not isinstance(basetype, File):
        return True
    return False



### PYTHONTOOL OUTPUTS ###
class PythonToolProcessOutputFactory:
    def __init__(self, out: TOutput, tool: PythonTool, sources: dict[str, Any]) -> None:
        self.out = out
        self.tool = tool
        self.sources = sources
        self.otype = get_otype(self.out)
        self.strategy_map = {
            OType.FILE: self.file_output,
            OType.NON_FILE: self.nonfile_output,
            OType.NON_FILE_ARRAY: self.nonfile_array_output,
        }

    @property
    def optional(self) -> bool:
        if self.out.outtype.optional:
            return True
        return False
    
    @property
    def target_file(self) -> str:
        return f'{settings.PYTHON_CODE_OUTPUT_FILENAME_PREFIX}{self.out.id()}'

    def create(self) -> ProcessOutput:
        strategy = self.strategy_map[self.otype]
        process_output = strategy()
        return process_output
    
    def file_output(self) -> PathProcessOutput:
        expr = f'"{self.target_file}"'
        new_output = PathProcessOutput(
            name=self.out.id(), 
            is_optional=self.optional, 
            expression=expr
        )
        return new_output
    
    def nonfile_output(self) -> ValProcessOutput:
        expr = f'\"${{file(\"${{task.workDir}}/{self.target_file}\").text}}"'
        new_output = ValProcessOutput(
            name=self.out.id(), 
            is_optional=self.optional, 
            expression=expr
        )
        return new_output

    def nonfile_array_output(self) -> ValProcessOutput:
        expr = f'\"${{file(\"${{task.workDir}}/{self.target_file}\").text.split(\',\')}}"'
        new_output = ValProcessOutput(
            name=self.out.id(), 
            is_optional=self.optional, 
            expression=expr
        )
        return new_output

    # def create_path_output_secondaries(self, ext: str) -> PathProcessOutput:
    #     # array of secondaries
    #     outputs: list[ProcessOutput] = []
    #     assert(isinstance(self.basetype, File))
    #     exts = secondaries.get_names(self.basetype)
    #     for ext in exts:
    #         outputs.append(self.create_path_output_secondaries(ext))
    #     return outputs
