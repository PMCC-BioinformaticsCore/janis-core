

from typing import Optional
from abc import ABC, abstractmethod
from dataclasses import dataclass

from .. import utils

from janis_core import ToolInput
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




def create_inputs(inp: ToolInput) -> list[ProcessInput]:
    datatype: DataType = inp.input_type
    if isinstance(datatype, Array):
        return create_inputs_array(inp)
    else:
        return create_inputs_single(inp)


def create_inputs_array(inp: ToolInput) -> list[ProcessInput]:
    basetype: Optional[DataType] = utils.get_base_type(inp.input_type)
    assert(basetype)

    # secondaries array
    if isinstance(basetype, File) and basetype.has_secondary_files():
        # a path input per file type
        inputs: list[ProcessInput] = []
        # get all extensions 
        exts = utils.get_extensions(basetype)
        for ext in exts:
            inputs.append(create_path_input_secondaries(inp, ext))
        return inputs

    # file array
    if isinstance(basetype, (File, Directory)):
        return [create_path_input(inp)]

    # nonfile array
    return [create_val_input(inp)]


def create_inputs_single(inp: ToolInput) -> list[ProcessInput]:
    basetype: Optional[DataType] = utils.get_base_type(inp.input_type)
    assert(basetype)
        
    # file secondaries
    if isinstance(basetype, File) and basetype.has_secondary_files():
        inputs = [create_tuple_input_secondaries(inp)]
        return inputs
    
    # file
    if isinstance(basetype, (File, Directory)):
        return [create_path_input(inp)]
    
    # nonfile
    return [create_val_input(inp)]


def create_path_input(inp: ToolInput) -> PathProcessInput:
    new_input = PathProcessInput(name=inp.id())
    new_input.presents_as = inp.presents_as
    return new_input

def create_val_input(inp: ToolInput) -> ValProcessInput:
    new_input = ValProcessInput(name=inp.id())
    return new_input

def create_path_input_secondaries(inp: ToolInput, ext: str) -> PathProcessInput:
    # TODO ignoring secondaries_presents_as for now!
    name = f'{inp.id()}_{ext}s'
    new_input = PathProcessInput(name=name)
    return new_input

def create_tuple_input_secondaries(inp: ToolInput) -> TupleProcessInput:
    assert(isinstance(inp.input_type, File))
    qualifiers: list[str] = []
    subnames: list[str] = []

    # tuple sub-element for each file
    exts = utils.get_extensions(inp.input_type)
    for ext in exts:
        qualifiers.append('path')
        subnames.append(ext)
    
    new_input = TupleProcessInput(
        name=inp.id(), 
        qualifiers=qualifiers, 
        subnames=subnames
    )
    return new_input