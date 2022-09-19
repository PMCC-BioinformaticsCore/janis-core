

from __future__ import annotations
from abc import ABC, abstractmethod
from dataclasses import dataclass
from io import RawIOBase
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
        width_col_1 = max([p.width for p in self.params])
        outstr = ''
        for p in self.params:
            outstr += f'{p.param:<{width_col_1}} = {p.default}\n'
        return outstr



class ParamDeclaration(NFBase):
    def __init__(self, name: str, default: Any):
        self.name = name
        self._default = default

    @property
    def default(self) -> str:
        if isinstance(self._default, list):
            default = [str(x) for x in self._default] # type: ignore
            default = ', '.join(default)  
        else:
            default = str(self._default)
        if default in formatting.type_keyword_map:
            default = formatting.type_keyword_map[default]
        return default
    
    @property
    def param(self) -> str:
        return f'params.{self.name}'

    @property
    def width(self) -> int:
        return len(self.param)
    
    def get_string(self) -> str:
        # leave this unimplemented! 
        # bad architecture means the method has to exist. 
        raise NotImplementedError  




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
