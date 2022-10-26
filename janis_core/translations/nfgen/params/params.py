

from __future__ import annotations
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any, Optional

from ..common import NFBase
from .. import utils



### ordering

@dataclass
class OrderingMethod(ABC):
    @abstractmethod
    def order(self, params: list[Param]) -> list[Param]:
        ...

@dataclass
class NotNullPriority(OrderingMethod):
    def order(self, params: list[Param]) -> list[Param]:
        params.sort(key=lambda x: x.value != 'null') 
        return params

@dataclass
class Alphabetical(OrderingMethod):
    def order(self, params: list[Param]) -> list[Param]:
        params.sort(key=lambda x: x.name) 
        return params

@dataclass
class MandatoryPriority(OrderingMethod):
    def order(self, params: list[Param]) -> list[Param]:
        params.sort(key=lambda x: x.optional or False) 
        return params

@dataclass
class FileTypePriority(OrderingMethod):
    def order(self, params: list[Param]) -> list[Param]:
        params.sort(key=lambda x: x.dtype == 'File', reverse=True) 
        return params

orderers: list[OrderingMethod] = [
    Alphabetical(),
    #NotNullPriority(),
    FileTypePriority(),
    MandatoryPriority()
]

### main classes 

class ParamRegister(NFBase):
    def __init__(self):
        self.params: dict[str, Param] = {}

    @property
    def ordered_params(self) -> list[Param]:
        params = list(self.params.values())
        for orderer in orderers:
            params = orderer.order(params)
        return params

    def get_string(self) -> str:
        width_col_1 = max([p.width for p in self.params.values()])
        outstr = ''
        for p in self.ordered_params:
            outstr += f'{p.text:<{width_col_1}} = {p.default}\n'
        return outstr


@dataclass
class Param(NFBase):
    varname: str
    scope: Optional[list[str]]=None
    dtype: Optional[str]=None
    optional: Optional[bool]=None
    default: Any=None

    @property
    def name(self) -> str:
        if self.scope and len(self.scope) > 0:
            name = f"{'_'.join(self.scope)}_{self.varname}"
        else:
            name = self.varname 
        return name.lower()

    @property
    def value(self) -> str:
        # get the default value as string
        if self.default == '':
            value = 'None'
        elif isinstance(self.default, list):
            value = [str(x) for x in self.default] # type: ignore
            value = ', '.join(value)  
        else:
            value = str(self.default)
        # cast 'None' to 'null' etc
        if value in utils.type_keyword_map:
            value = utils.type_keyword_map[value]
        return value
    
    @property
    def text(self) -> str:
        return f'params.{self.name}'
    
    @property
    def width(self) -> int:
        return len(self.text)
    
    def get_string(self) -> str:
        # leave this unimplemented! 
        # architecture mandates the method has to exist, but not needed. 
        raise NotImplementedError  



### module entry points

register = ParamRegister()

def add(varname: str, 
        scope: Optional[list[str]]=None, 
        dtype: Optional[str]=None, 
        optional: Optional[bool]=None, 
        default: Any=None
    ) -> None:
    global register
    if exists(varname, scope):
        raise RuntimeError(f'param already exists. scope: {scope}, varname: {varname}')
    else:
        param = Param(varname, scope, dtype, optional, default)
        register.params[param.name] = param

def update(varname: str, 
        scope: Optional[list[str]]=None, 
        dtype: Optional[str]=None, 
        optional: Optional[bool]=None, 
        default: Any=None
    ) -> None:
    global register
    param = Param(varname, scope, dtype, optional, default)
    register.params[param.name] = param
    
def exists(varname: str, scope: Optional[list[str]]=None) -> bool:
    global register
    param = Param(varname, scope)
    if param.name in register.params:
        return True
    return False

def get(varname: str, scope: Optional[list[str]]=None) -> Param:
    global register
    param = Param(varname, scope)
    return register.params[param.name]

def in_scope(scope: list[str]) -> list[Param]:
    global register
    return [x for x in register.ordered_params if x.scope == scope]

def getall() -> list[Param]:
    global register
    return register.ordered_params

def getstr() -> str:
    global register
    return register.get_string()

def serialize() -> dict[str, Any]:
    global register
    the_dict: dict[str, Any] = {}
    for p in getall():
        the_dict[p.name] = p.value
    return the_dict





