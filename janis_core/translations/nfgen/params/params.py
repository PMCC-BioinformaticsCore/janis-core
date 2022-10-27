

from __future__ import annotations
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any, Optional

from janis_core.workflow.workflow import InputNode
from janis_core.workflow.workflow import StepNode
from janis_core.workflow.workflow import Workflow
from janis_core.tool.commandtool import CommandTool
from janis_core.tool.commandtool import ToolInput

from ..common import NFBase
from .. import utils


def register(
    the_entity: Optional[Workflow | StepNode | CommandTool]=None,
    the_dict: Optional[dict[str, Any]]=None, 
    scope: Optional[list[str]]=None, 
    override: bool=True
    ) -> None:
    if not isinstance(scope, list):
        scope = []
    if the_entity:
        register_params_for_entity(the_entity, scope, override)
    elif the_dict:
        register_params_for_dict(the_dict, scope, override)
    else:
        raise RuntimeError("nothing to register: please supply 'the_entity' or 'the_dict'")

def register_params_for_entity(
    entity: Workflow | StepNode | CommandTool, 
    scope: list[str], 
    override: bool
    ) -> None:
    """
    Registers a param tool or workflow inputs.
    Workflow / tool inputs which are exposed to the user must to be listed
    as part of the global params object. 
    """
    if isinstance(entity, Workflow):
        inputs = utils.get_workflow_inputs(entity)
        for inp in inputs:
            register_wfinp_param(inp, scope, override=override)
    elif isinstance(entity, StepNode):
        inputs = utils.get_exposed_tool_inputs(entity.tool, entity.sources)
        for inp in inputs:
            src = entity.sources[inp.id()]
            val = utils.get_source_value(src)
            register_toolinp_param(inp, scope=scope, value=val, override=override)
        print('\n\nPARAMS ----')
        print(getstr())
    elif isinstance(entity, CommandTool):  # type: ignore
        inputs = utils.get_exposed_tool_inputs(entity, {})
        for inp in inputs:
            register_toolinp_param(inp, scope=scope, override=override)
    else:
        raise RuntimeError


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
            dtype=inp.datatype.name(),
            optional=inp.datatype.optional,
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
            dtype=inp.input_type.name(),
            optional=inp.input_type.optional,
            scope=scope,
            default=utils.wrap_value(default, inp),
            is_wf_input=False
        )

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

@dataclass
class WfInputPriority(OrderingMethod):
    def order(self, params: list[Param]) -> list[Param]:
        params.sort(key=lambda x: x.is_wf_input, reverse=True) 
        return params

orderers: list[OrderingMethod] = [
    #NotNullPriority(),
    #FileTypePriority(),
    #MandatoryPriority(),
    Alphabetical(),
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
    is_wf_input: bool=False

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

param_register = ParamRegister()

def add(varname: str, 
        scope: Optional[list[str]]=None, 
        dtype: Optional[str]=None, 
        optional: Optional[bool]=None, 
        default: Any=None,
        is_wf_input: bool=False
        
    ) -> None:
    global param_register
    param = Param(varname, scope, dtype, optional, default, is_wf_input)
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

def getstr() -> str:
    global param_register
    return param_register.get_string()

def serialize() -> dict[str, Any]:
    global param_register
    the_dict: dict[str, Any] = {}
    for p in getall():
        the_dict[p.name] = p.value
    return the_dict





