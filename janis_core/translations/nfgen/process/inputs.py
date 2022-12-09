

from typing import Optional
from abc import ABC, abstractmethod
from dataclasses import dataclass

from .. import nfgen_utils
from .. import secondaries

from janis_core import ToolInput, TInput
from janis_core.types import File, Directory, Array, DataType


@dataclass
class ProcessInput(ABC):
    name: str
    
    @abstractmethod
    def get_string(self) -> str:
        ...


@dataclass
class ValProcessInput(ProcessInput):

    def get_string(self) -> str:
        return f'val {self.name}'


@dataclass
class PathProcessInput(ProcessInput):
    presents_as: Optional[str]=None

    @property
    def stage_as(self) -> str:
        if self.presents_as:
            return f", stageAs: '{self.presents_as}'"
        return ''
    
    def get_string(self) -> str:
        return f'path {self.name}{self.stage_as}'


@dataclass
class TupleProcessInput(ProcessInput):
    qualifiers: list[str]
    subnames: list[str]

    @property
    def fields(self) -> str:
        out: str = ''
        for qual, subname in zip(self.qualifiers, self.subnames):
            out += f'{qual}({subname}), '
        out = out.rstrip(', ') # strip off the last comma & space
        return out

    def get_string(self) -> str:
        return f'tuple {self.fields}'




def create_inputs(inp: ToolInput | TInput) -> list[ProcessInput]:
    dtype: DataType = inp.input_type if isinstance(inp, ToolInput) else inp.intype # type: ignore
    datatype: DataType = dtype
    if isinstance(datatype, Array):
        return create_inputs_array(inp)
    else:
        return create_inputs_single(inp)

def create_inputs_array(inp: ToolInput | TInput) -> list[ProcessInput]:
    dtype: DataType = inp.input_type if isinstance(inp, ToolInput) else inp.intype # type: ignore
    basetype: Optional[DataType] = nfgen_utils.get_base_type(dtype)
    assert(basetype)

    # secondaries array
    if isinstance(basetype, File) and basetype.has_secondary_files():
        # a path input per file type
        inputs: list[ProcessInput] = []
        # get all extensions 
        exts = secondaries.get_names(basetype)
        for ext in exts:
            inputs.append(create_path_input_secondaries(inp, ext))
        return inputs

    # file array
    if isinstance(basetype, (File, Directory)):
        return [create_path_input(inp)]

    # nonfile array
    return [create_val_input(inp)]


def create_inputs_single(inp: ToolInput | TInput) -> list[ProcessInput]:
    dtype: DataType = inp.input_type if isinstance(inp, ToolInput) else inp.intype # type: ignore
    basetype: Optional[DataType] = nfgen_utils.get_base_type(dtype)
    assert(basetype)
        
    # file secondaries
    if isinstance(basetype, File) and basetype.has_secondary_files():
        inputs = [create_tuple_input_secondaries(inp)]
        return inputs # type: ignore
    
    # file
    if isinstance(basetype, (File, Directory)):
        return [create_path_input(inp)]
    
    # nonfile
    return [create_val_input(inp)]


def create_path_input(inp: ToolInput | TInput) -> PathProcessInput:
    new_input = PathProcessInput(name=inp.id())
    new_input.presents_as = None
    if isinstance(inp, ToolInput):
        new_input.presents_as = inp.presents_as
    return new_input

def create_val_input(inp: ToolInput | TInput) -> ValProcessInput:
    new_input = ValProcessInput(name=inp.id())
    return new_input

def create_path_input_secondaries(inp: ToolInput | TInput, ext: str) -> PathProcessInput:
    # TODO ignoring secondaries_presents_as for now!
    name = f'{inp.id()}_{ext}s'
    new_input = PathProcessInput(name=name)
    return new_input

def create_tuple_input_secondaries(inp: ToolInput | TInput) -> TupleProcessInput:
    dtype: DataType = inp.input_type if isinstance(inp, ToolInput) else inp.intype # type: ignore
    assert(isinstance(dtype, File))
    qualifiers: list[str] = []
    subnames: list[str] = []

    # tuple sub-element for each file
    exts = secondaries.get_names(dtype)
    for ext in exts:
        qualifiers.append('path')
        subnames.append(ext)
    
    new_input = TupleProcessInput(
        name=inp.id(), 
        qualifiers=qualifiers, 
        subnames=subnames
    )
    return new_input