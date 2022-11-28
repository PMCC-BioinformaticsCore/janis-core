

from __future__ import annotations
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any, Optional

from janis_core.types import (
    DataType,
    String, 
    Array, 
    File
)

from janis_core.translations.nfgen import NFBase
from janis_core.translations.nfgen import nfgen_utils
from janis_core.translations.nfgen.casefmt import to_case
from janis_core.translations.nfgen import settings


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
            basetype = nfgen_utils.get_base_type(p.dtype)
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


class ParamRegister(NFBase):
    def __init__(self):
        self.params: list[Param] = []

    @property
    def ordered_params(self) -> list[Param]:
        params = self.params
        for orderer in orderers:
            params = orderer.order(params)
        return params
    
    def get_string(self) -> str:
        # leave this unimplemented! 
        # architecture mandates the method has to exist, but not needed. 
        raise NotImplementedError  


@dataclass
class Param(NFBase):
    var_name: str
    var_scope: list[str]
    dtype: Optional[DataType]=None
    default: Any=None
    is_channel_input: bool=False
    name_override: Optional[str]=None
    janis_uuid: Optional[str]=None

    # def __post_init__(self):
    #     self.uuid = str(uuid4())

    @property
    def name(self) -> str:
        if self.name_override:
            base = to_case(self.name_override, settings.NEXTFLOW_PARAM_CASE)
        else:
            base = to_case(self.var_name, settings.NEXTFLOW_PARAM_CASE)
        if self.var_scope:
            scope = [to_case(x, settings.NEXTFLOW_PARAM_CASE) for x in self.var_scope]
            name = f"{'.'.join(scope)}.{base}"
        else:
            name = base
        return name

    @property
    def groovy_value(self) -> str:
        # get the default value as groovy code string
        # TODO I am dubious about this
        # if self.default == '':
        #     val = None
        if isinstance(self.dtype, Array) and self.default is None:
            val: list[str] = []
        else:
            val = self.default
        return nfgen_utils.to_groovy(val, self.dtype)
    
    @property
    def width(self) -> int:
        return len(self.name)

    def get_string(self) -> str:
        # leave this unimplemented! 
        # architecture mandates the method has to exist, but not needed. 
        raise NotImplementedError  


### instantiation of param register & default params

param_register = ParamRegister()

default_params = [
    Param(
        var_name='outdir',
        var_scope=[],
        dtype=String(),
        default='"outputs"',
        is_channel_input=False,
    )
]

for param in default_params:
    param_register.params.append(param)


### MODULE ENTRY POINTS

def add(var_name: str, 
        var_scope: Optional[list[str]]=None, 
        dtype: Optional[DataType]=None, 
        default: Any=None,
        is_channel_input: bool=False,
        name_override: Optional[str]=None,
        janis_uuid: Optional[str]=None
    ) -> Param:
    global param_register
    var_scope = var_scope if var_scope else []
    param = Param(
        var_name=var_name, 
        var_scope=var_scope, 
        dtype=dtype, 
        default=default, 
        is_channel_input=is_channel_input, 
        name_override=name_override,
        janis_uuid=janis_uuid
    )
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


# def add_link(upstream_janis_uuid: str, downstream_janis_uuid: str) -> None:
#     global param_register
#     assert(exists(upstream_janis_uuid))
#     param_register.upstream[downstream_janis_uuid].append(upstream_janis_uuid)
#     param_register.downstream[upstream_janis_uuid].append(downstream_janis_uuid)

# def upstream(query: str) -> Optional[Param]:
#     if query in param_register.upstream:
#         upstream_janis_uuid = param_register.upstream[query][0]
#         return get(upstream_janis_uuid)
#     return None

# def upstream_multi(query: str) -> list[Param]:
#     params: list[Param] = []
#     for upstream_janis_uuid in param_register.upstream[query]:
#         params += getall(upstream_janis_uuid)
#     return params








