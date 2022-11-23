

from __future__ import annotations
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any, Optional

from janis_core.workflow.workflow import InputNode
from janis_core.workflow.workflow import Workflow
from janis_core.tool.commandtool import CommandTool
from janis_core.types import (
    DataType,
    String, 
    Array, 
    File
)

from janis_core.translations.nfgen import NFBase
from janis_core.translations.nfgen import utils
from janis_core.translations.nfgen.casefmt import to_case
from janis_core.translations.nfgen import settings
from uuid import uuid4




"""
MINIMAL PROCESS
- wf inputs: param for all wf inputs
- tool inputs: param for non-process-inputs fed value using step.sources

FULL PROCESS
- wf inputs: param for all wf inputs
- tool inputs: param for all non-process-inputs
"""




# # tool inputs
# def register_toolinp_param(
#     inp: ToolInput, 
#     scope: Optional[list[str]],
#     default: Optional[Any]='__UwU_PlaceholdeR_UwU__',
#     override: bool=False
#     ) -> None:
#     if (not exists(inp.id(), scope)) or (exists(inp.id(), scope) and override):
#         # valid 'value' includes 'None', so must use placeholder
#         if default != '__UwU_PlaceholdeR_UwU__':
#             default = default
#         else:
#             default = inp.default 
#         add(
#             varname=inp.id(),
#             reference=inp.id(),
#             dtype=inp.input_type,
#             scope=scope,
#             default=default,
#             is_wf_input=False
#         )



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
            dtype = p.dtype
            while isinstance(dtype, Array):
                dtype = dtype.subtype()
            if isinstance(dtype, File):
                top.append(p)
            else:
                bottom.append(p)
        return top + bottom

@dataclass
class WfInputPriority(OrderingMethod):
    def order(self, params: list[Param]) -> list[Param]:
        return sorted(params, key=lambda x: x.is_wf_input, reverse=True) 

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
        self.params: dict[str, Param] = {}

    @property
    def ordered_params(self) -> list[Param]:
        params = list(self.params.values())
        for orderer in orderers:
            params = orderer.order(params)
        return params
    
    def get_string(self) -> str:
        # leave this unimplemented! 
        # architecture mandates the method has to exist, but not needed. 
        raise NotImplementedError  


@dataclass
class Param(NFBase):
    ref_name: str
    ref_scope: list[str]
    dtype: Optional[DataType]=None
    default: Any=None
    is_wf_input: bool=False
    name_override: Optional[str]=None

    def __post_init__(self):
        self.uuid = uuid4()

    @property
    def name(self) -> str:
        if self.name_override:
            base = self.name_override
        else:
            base = self.ref_name
        if self.ref_scope:
            full = f"{'_'.join(self.ref_scope)}_{base}"
        else:
            full = base
        name = to_case(full, settings.NEXTFLOW_PARAM_CASE)
        return name

    @property
    def groovy_value(self) -> str:
        # get the default value as string
        # TODO I am dubious about this
        # if self.default == '':
        #     val = None
        if isinstance(self.dtype, Array) and self.default is None:
            val: list[str] = []
        else:
            val = self.default
        return utils.to_groovy(val, self.dtype)
    
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
        ref_name='outdir',
        ref_scope=[],
        dtype=String(),
        default='"outputs"',
        is_wf_input=False,
    )
]

for param in default_params:
    param_register.params[param.name] = param


### MODULE ENTRY POINTS

def add(ref_name: str, 
        ref_scope: Optional[list[str]]=None, 
        dtype: Optional[DataType]=None, 
        default: Any=None,
        is_wf_input: bool=False,
        name_override: Optional[str]=None
    ) -> Param:
    global param_register
    ref_scope = ref_scope if ref_scope else []
    param = Param(
        ref_name=ref_name, 
        ref_scope=ref_scope, 
        dtype=dtype, 
        default=default, 
        is_wf_input=is_wf_input, 
        name_override=name_override
    )
    param_register.params[param.name] = param
    return param
    
def exists(name: str, scope: Optional[list[str]]=None) -> bool:
    global param_register
    scope = scope if scope else []
    for param in param_register.ordered_params:
        if param.ref_name == name and param.ref_scope == scope:
            return True
    return False

def get(name: str, scope: Optional[list[str]]=None) -> Param:
    global param_register
    return getall(name, scope)[0]

def getall(name: Optional[str]=None, scope: Optional[list[str]]=None) -> list[Param]:
    global param_register
    params = param_register.ordered_params
    if name:
        params = [x for x in params if x.ref_name == name]
    if scope:
        params = [x for x in params if x.ref_scope == scope]
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






