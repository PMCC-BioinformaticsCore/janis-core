

from __future__ import annotations
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any

from .common import NFBase
from . import utils



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
    def __init__(self, name: str, default: Any=None):
        self.name = name
        self._default = default

    @property
    def default(self) -> str:
        if isinstance(self._default, list):
            default = [str(x) for x in self._default] # type: ignore
            default = ', '.join(default)  
        else:
            default = str(self._default)
        if default in utils.type_keyword_map:
            default = utils.type_keyword_map[default]
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
