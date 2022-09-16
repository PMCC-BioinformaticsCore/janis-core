

from __future__ import annotations
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any
from .common import NFBase
from . import formatting


@dataclass
class OrderingMethod(ABC):
    @abstractmethod
    def order(self, params: list[ParamDeclaration]) -> list[ParamDeclaration]:
        ...

@dataclass
class NotNullPriority(OrderingMethod):
    def order(self, params: list[ParamDeclaration]) -> list[ParamDeclaration]:
        params.sort(key=lambda x: x.default != 'null', reverse=True) 
        return params

@dataclass
class Alphabetical(OrderingMethod):
    def order(self, params: list[ParamDeclaration]) -> list[ParamDeclaration]:
        params.sort(key=lambda x: x.name) 
        return params

orderers: list[OrderingMethod] = [
    Alphabetical(),
    NotNullPriority(),
]



class ParamDeclarationBlock(NFBase):
    def __init__(self, params: list[ParamDeclaration]):
        self.params = params

    @property
    def ordered_params(self) -> list[ParamDeclaration]:
        params = self.params
        for orderer in orderers:
            params = orderer.order(params)
        return params

    def get_string(self) -> str:
        return '\n'.join([ch.get_string() for ch in self.ordered_params])
    



class ParamDeclaration(NFBase):
    def __init__(self, name: str, default: Any):
        self.name = name
        self._default = default

    @property
    def default(self) -> str:
        default = str(self._default)
        if default in formatting.type_keyword_map:
            default = formatting.type_keyword_map[default]
        return default
    
    def get_string(self) -> str:
        if isinstance(self.default, list):
            return f'params.{self.name} = {", ".join(self.default)}'  
        else:
            return f'params.{self.name} = {self.default}'



    # # has default
    # if task_input.default is not None:
    #     if method == 'value' and utils.is_path(task_input):
    #         sources = [f"file('{task_input.default}')"]
    #     elif isinstance(task_input.default, list):
    #         sources = ', '.join(task_input.default)
    #     else:
    #         sources = [task_input.default]
    # else:
    #     sources = ['None']
    # return sources
