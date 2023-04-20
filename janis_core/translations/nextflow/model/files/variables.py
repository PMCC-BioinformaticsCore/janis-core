

from __future__ import annotations
from dataclasses import dataclass

from abc import ABC, abstractmethod

from janis_core.types import DataType, File, Filename, Directory
from janis_core import translation_utils as utils


### main class 

@dataclass
class NFVariableDefinition:
    name: str
    source: str
    dtype: DataType

    @property
    def name_width(self) -> int:
        return len(self.name)

    def get_string(self) -> str:
        return f'file( {self.source} )'
    

@dataclass
class NFVariableDefinitionBlock:
    variables: list[NFVariableDefinition]

    @property
    def def_width(self) -> int:
        return max(v.name_width for v in self.variables) + 1

    def get_string(self) -> str:
        out: str = ''
        out += '// data which will be passed as optional files\n'
        ordered_variables = order_nf_variables(self.variables)
        for v in ordered_variables:
            line = f'{v.name:<{self.def_width}} = {v.get_string()}\n'
            out += line
        return out


### ordering 

@dataclass
class OrderingMethod(ABC):
    @abstractmethod
    def order(self, variables: list[NFVariableDefinition]) -> list[NFVariableDefinition]:
        ...

@dataclass
class FileTypePriority(OrderingMethod):
    def order(self, variables: list[NFVariableDefinition]) -> list[NFVariableDefinition]:
        top: list[NFVariableDefinition] = []
        bottom: list[NFVariableDefinition] = []
        for p in variables:
            basetype = utils.get_base_type(p.dtype)
            basetype = utils.ensure_single_type(basetype)
            if isinstance(basetype, (File, Filename, Directory)):
                top.append(p)
            else:
                bottom.append(p)
        return top + bottom

@dataclass
class Alphabetical(OrderingMethod):
    def order(self, variables: list[NFVariableDefinition]) -> list[NFVariableDefinition]:
        return sorted(variables, key=lambda x: x.name) 

orderers: list[OrderingMethod] = [
    FileTypePriority(),
    Alphabetical(),
]

def order_nf_variables(variables: list[NFVariableDefinition]) -> list[NFVariableDefinition]:
    for orderer in orderers:
        variables = orderer.order(variables)
    return variables


