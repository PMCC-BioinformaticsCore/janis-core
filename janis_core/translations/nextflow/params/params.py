

from __future__ import annotations
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any, Optional

from janis_core import settings
from janis_core.types import (
    DataType,
    Directory,
    Array, 
    File
)

from janis_core import translation_utils as utils
from .. import nfgen_utils
from .. import naming
from ..scope import Scope


### ORDERING

@dataclass
class OrderingMethod(ABC):
    @abstractmethod
    def order(self, params: list[Param]) -> list[Param]:
        ...

@dataclass
class NotNullPriority(OrderingMethod):
    def order(self, params: list[Param]) -> list[Param]:
        return sorted(params, key=lambda x: x.groovy_value != 'null') 

@dataclass
class Alphabetical(OrderingMethod):
    def order(self, params: list[Param]) -> list[Param]:
        return sorted(params, key=lambda x: x.name) 

@dataclass
class MandatoryPriority(OrderingMethod):
    def order(self, params: list[Param]) -> list[Param]:
        top: list[Param] = []
        bottom: list[Param] = []
        for p in params:
            if p.janis_type and p.janis_type.optional == False:
                top.append(p)
            else:
                bottom.append(p)
        return top + bottom

@dataclass
class FileTypePriority(OrderingMethod):
    def order(self, params: list[Param]) -> list[Param]:
        top: list[Param] = []
        bottom: list[Param] = []
        for p in params:
            basetype = utils.get_base_type(p.janis_type)
            basetype = utils.ensure_single_type(basetype)
            if isinstance(basetype, File):
                top.append(p)
            else:
                bottom.append(p)
        return top + bottom

@dataclass
class WfInputPriority(OrderingMethod):
    def order(self, params: list[Param]) -> list[Param]:
        return sorted(params, key=lambda x: x.is_channel_input, reverse=True) 

orderers: list[OrderingMethod] = [
    #NotNullPriority(),
    #FileTypePriority(),
    #MandatoryPriority(),
    #Alphabetical(),
    WfInputPriority(),
]


### MAIN CLASSES  


class ParamRegister:
    def __init__(self):
        self.params: list[Param] = []

    @property
    def ordered_params(self) -> list[Param]:
        params = self.params
        for orderer in orderers:
            params = orderer.order(params)
        return params
    
    def get_string(self) -> str:
        raise NotImplementedError  


@dataclass
class Param:
    name: str
    scope: Scope
    default: Any=None
    is_channel_input: bool=False
    janis_type: Optional[DataType]=None
    janis_uuid: Optional[str]=None

    @property
    def groovy_value(self) -> str:
        # get the default value as groovy code string
        if isinstance(self.janis_type, Array) and self.default is None:
            val: list[str] = []
        else:
            val = self.default
        return nfgen_utils.to_groovy(val, self.janis_type)
    
    @property
    def width(self) -> int:
        return len(self.name)

    def get_string(self) -> str:
        raise NotImplementedError  



### MODULE ENTRY POINTS

def add(janis_tag: str, 
        scope: Scope,
        default: Any=None,
        is_channel_input: bool=False,
        name_override: Optional[str]=None,
        janis_dtype: Optional[DataType]=None, 
        janis_uuid: Optional[str]=None
    ) -> Param:
    global param_register
    # channel name
    name = naming.constructs.gen_varname_param(janis_tag, scope, name_override, janis_dtype)
    # create param
    param = Param(name, scope, default, is_channel_input, janis_dtype, janis_uuid)
    # register param
    param_register.params.append(param)
    return param
    
def exists(janis_uuid: str) -> bool:
    global param_register
    for param in param_register.ordered_params:
        if param.janis_uuid == janis_uuid:
            return True
    return False

def get(janis_uuid: str) -> Param:
    return getall(janis_uuid)[0]

def getall(janis_uuid: Optional[str]=None) -> list[Param]:
    global param_register
    params = param_register.ordered_params
    if janis_uuid:
        params = [x for x in params if x.janis_uuid == janis_uuid]
    return params

def serialize() -> dict[str, Any]:
    global param_register
    the_dict: dict[str, Any] = {}
    for p in getall():
        the_dict[p.name] = p.groovy_value
    return the_dict

def clear() -> None:
    global param_register 
    param_register = ParamRegister()
    add_default_params()



### instantiation of param register & default params

default_params = [
    {
        'janis_tag': None,
        'scope': Scope(),
        'default': './outputs',
        'is_channel_input': False,
        'name_override': 'outdir',
        'janis_dtype': Directory(),
        'janis_uuid': None
    }
]

def add_default_params():
    for p in default_params:
        add(
            janis_tag=p['janis_tag'], 
            scope=p['scope'],
            default=p['default'],
            is_channel_input=p['is_channel_input'],
            name_override=p['name_override'],
            janis_dtype=p['janis_dtype'],
            janis_uuid=p['janis_uuid'],
        )

param_register = ParamRegister()
add_default_params()