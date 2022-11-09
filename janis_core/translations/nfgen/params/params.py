

from __future__ import annotations
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any, Optional

from janis_core.workflow.workflow import InputNode
from janis_core.workflow.workflow import Workflow
from janis_core.tool.commandtool import CommandTool
from janis_core.tool.commandtool import ToolInput
from janis_core.types import (
    DataType,
    String, 
    Array, 
    File
)

from ..common import NFBase
from .. import utils



"""
MINIMAL PROCESS
- wf inputs: param for all wf inputs
- tool inputs: param for non-process-inputs fed value using step.sources

FULL PROCESS
- wf inputs: param for all wf inputs
- tool inputs: param for all non-process-inputs
"""


def register(
    the_entity: Optional[Workflow | CommandTool]=None,
    the_dict: Optional[dict[str, Any]]=None, 
    sources: Optional[dict[str, Any]]=None, 
    scope: Optional[list[str]]=None, 
    override: bool=True
    ) -> None:

    if not isinstance(scope, list):
        scope = []
    if isinstance(the_entity, Workflow):
        register_params_for_wf_inputs(the_entity, scope, override, sources)
    elif isinstance(the_entity, CommandTool):
        register_params_for_tool(the_entity, scope, override, sources)
    elif the_dict is not None:
        register_params_for_dict(the_dict, scope, override)
    else:
        raise RuntimeError("nothing to register: please supply 'the_entity' or 'the_dict'")


def register_params_for_wf_inputs(
    workflow: Workflow,
    scope: list[str], 
    override: bool,
    sources: Optional[dict[str, Any]]=None, 
    ) -> None:
    param_ids = utils.get_channel_input_ids(workflow)
    param_inputs = utils.items_with_id(list(workflow.input_nodes.values()), param_ids)
    for inp in param_inputs:
        val = None
        if sources is not None and inp.id() in sources:
            src = sources[inp.id()]
            val = utils.get_source_value(src)
        register_wfinp_param(inp, scope=scope, value=val, override=override)

def register_params_for_tool(
    tool: CommandTool, 
    scope: list[str], 
    override: bool,
    sources: Optional[dict[str, Any]]=None, 
    ) -> None:
    """
    Registers a param tool or workflow inputs.
    Workflow / tool inputs which are exposed to the user must to be listed
    as part of the global params object. 
    """
    sources = sources if sources is not None else {}
    param_ids = utils.get_param_input_ids(tool, sources=sources)
    param_inputs = utils.items_with_id(tool.inputs(), param_ids)
    for inp in param_inputs:
        val = None
        if sources is not None and inp.id() in sources:
            src = sources[inp.id()]
            val = utils.get_source_value(src)
        register_toolinp_param(inp, scope=scope, value=val, override=override)

def register_params_for_dict(
    the_dict: dict[str, Any],
    scope: list[str],
    override: bool
    ) -> None:
    for k, v in the_dict.items():
        if (not exists(k, scope)) or (exists(k, scope) and override):
            add(
                varname=k,
                scope=scope,
                default=v,
                is_wf_input=False
            )


def register_wfinp_param(
    inp: InputNode, 
    scope: Optional[list[str]], 
    value: Optional[Any]=None, 
    override: bool=False
    ) -> None:
    varname = inp.id()
    # brackets just to be specific about the condition
    if (not exists(varname, scope)) or (exists(varname, scope) and override):
        default = value if value else utils.wrap_value(inp.default, inp)  # type: ignore
        add(
            varname=varname,
            scope=scope,
            dtype=inp.datatype,
            default=default,
            is_wf_input=True
        )

def register_toolinp_param(
    inp: ToolInput, 
    scope: Optional[list[str]],
    value: Optional[Any]='__UwU_PlaceholdeR_UwU__',
    override: bool=False
    ) -> None:
    varname = inp.id()
    if (not exists(varname, scope)) or (exists(varname, scope) and override):
        # valid 'value' includes 'None', so must use placeholder
        if value != '__UwU_PlaceholdeR_UwU__':
            default = value
        else:
            default = inp.default 
        add(
            varname=varname,
            dtype=inp.input_type,
            scope=scope,
            default=utils.wrap_value(default, inp),
            is_wf_input=False
        )



### ordering

@dataclass
class OrderingMethod(ABC):
    @abstractmethod
    def order(self, params: list[Param]) -> list[Param]:
        ...

@dataclass
class NotNullPriority(OrderingMethod):
    def order(self, params: list[Param]) -> list[Param]:
        return sorted(params, key=lambda x: x.value != 'null') 

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
        # leave this unimplemented! 
        # architecture mandates the method has to exist, but not needed. 
        raise NotImplementedError  


@dataclass
class Param(NFBase):
    varname: str
    scope: Optional[list[str]]=None
    dtype: Optional[DataType]=None
    default: Any=None
    is_wf_input: bool=False

    @property
    def name(self) -> str:
        if self.scope:
            name = f"{'_'.join(self.scope)}_{self.varname}"
        else:
            name = self.varname
        return name
        # return name.lower()

    @property
    def value(self) -> str:
        # get the default value as string
        # TODO I am dubious about this
        if self.default == '':
            value = 'None'
        if isinstance(self.default, list):
            value = [str(x) for x in self.default] # type: ignore
            value = ', '.join(value)  
        else:
            value = str(self.default)
        # cast 'None' to 'null' etc
        if value in utils.type_keyword_map:
            value = utils.type_keyword_map[value]
        return value
    
    @property
    def width(self) -> int:
        return len(self.name)
    
    def get_string(self) -> str:
        # leave this unimplemented! 
        # architecture mandates the method has to exist, but not needed. 
        raise NotImplementedError  


### instantiation of register & default params

param_register = ParamRegister()

default_params = [
    Param(
        varname='outdir',
        dtype=String(),
        default='"outputs"',
        is_wf_input=False,
    )
]

for param in default_params:
    param_register.params[param.name] = param


### module entry points

def add(varname: str, 
        scope: Optional[list[str]]=None, 
        dtype: Optional[DataType]=None, 
        default: Any=None,
        is_wf_input: bool=False
        
    ) -> None:
    global param_register
    param = Param(varname, scope, dtype, default, is_wf_input)
    param_register.params[param.name] = param
    
def exists(varname: str, scope: Optional[list[str]]=None) -> bool:
    global param_register
    param = Param(varname, scope)
    if param.name in param_register.params:
        return True
    return False

def get(varname: str, scope: Optional[list[str]]=None) -> Param:
    global param_register
    param = Param(varname, scope)
    return param_register.params[param.name]

def getall() -> list[Param]:
    global param_register
    return param_register.ordered_params

def in_scope(scope: list[str]) -> list[Param]:
    global param_register
    return [x for x in param_register.ordered_params if x.scope == scope]

def serialize() -> dict[str, Any]:
    global param_register
    the_dict: dict[str, Any] = {}
    for p in getall():
        the_dict[p.name] = p.value
    return the_dict

def clear() -> None:
    global param_register 
    param_register = ParamRegister()






