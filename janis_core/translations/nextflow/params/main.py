

from __future__ import annotations
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any, Optional

from janis_core.types import (
    DataType,
    Directory,
    Array, 
    File
)

from janis_core import translation_utils as utils
from .. import nfgen_utils
from .. import naming


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
            if p.dtype and p.dtype.optional == False:
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
            basetype = utils.get_base_type(p.dtype)
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
    #MandatoryPriority(),
    Alphabetical(),
    FileTypePriority(),
    # WfInputPriority(),
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
    
    def exists(self, tinput_id: str, task_id: str) -> bool:
        for p in self.params:
            if p.tinput_id == tinput_id and p.task_id == task_id:
                return True
        return False
    
    def get(self, tinput_id: str, task_id: str) -> Param:
        for p in self.params:
            if p.tinput_id == tinput_id and p.task_id == task_id:
                return p
        raise RuntimeError
    
    def get_string(self) -> str:
        raise NotImplementedError  


@dataclass
class Param:
    name: str
    tinput_id: str
    task_id: str
    default: Any
    dtype: Optional[DataType]

    @property
    def groovy_value(self) -> str:
        # get the default value as groovy code string
        if isinstance(self.dtype, Array) and self.default is None:
            val: list[str] = []
        else:
            val = self.default
        return nfgen_utils.to_groovy(val, self.dtype)
    
    @property
    def width(self) -> int:
        return len(self.name)

    def get_string(self) -> str:
        raise NotImplementedError  



### MODULE ENTRY POINTS

def add(
    task_id: str,
    tinput_id: str,
    name_override: Optional[str]=None,
    janis_dtype: Optional[DataType]=None,
    default: Optional[Any]=None,
    is_subtask_param: bool=False
) -> Param:
    global param_register
    name = naming.constructs.gen_varname_param(task_id, tinput_id, name_override, is_subtask_param)
    assert(tinput_id)
    param = Param(name, tinput_id, task_id, default, janis_dtype)
    param_register.params.append(param)
    return param

def exists(tinput_id: str, task_id: str) -> bool:
    global param_register
    return param_register.exists(tinput_id, task_id)

def get(tinput_id: str, task_id: str) -> Param:
    global param_register
    return param_register.get(tinput_id, task_id)

def getall() -> list[Param]:
    global param_register
    return param_register.ordered_params

def getstr() -> str:
    global param_register
    params = param_register.ordered_params
    param_names = [p.name for p in params]
    return '\n'.join(param_names)

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
        'task_id': '__DEFAULTS__',
        'tinput_id': '__DEFAULTS__',
        'name_override': 'outdir',
        'janis_dtype': Directory(),
        'default': './outputs',
        'is_subtask_param': False,
    }
]    

def add_default_params():
    for p in default_params:
        add(
            task_id=p['task_id'], 
            tinput_id=p['tinput_id'], 
            default=p['default'],
            name_override=p['name_override'],
            janis_dtype=p['janis_dtype'],
            is_subtask_param=p['is_subtask_param'],
        )

param_register = ParamRegister()
add_default_params()
