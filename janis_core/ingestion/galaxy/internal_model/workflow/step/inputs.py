
import ast 
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any, Optional

from janis_core.ingestion.galaxy.gxtool.command.components import OutputComponent
from janis_core.ingestion.galaxy.gxtool.command.components import InputComponent

from janis_core.ingestion.galaxy import tags
from janis_core.ingestion.galaxy import expressions


@dataclass
class InputValue(ABC):
    component: InputComponent

    def __post_init__(self):
        self.scatter: bool = False

    @property
    def comptype(self) -> str:
        return type(self.component).__name__.lower() 

    @property
    def input_tag(self) -> str:
        """get the str tag for this tool input"""
        return self.component.tag
    
    @property
    @abstractmethod
    def wrapped_value(self) -> str:
        """
        get the str value for this tool input.
        for workflow inputs & connections, returns janis style
        workflow object name references.
        """
        ...


def get_comptype(component: InputComponent | OutputComponent) -> str:
    return type(component).__name__.lower() 

def is_bool(value: str) -> bool:
    if isinstance(value, bool):
        return True
    return False

def is_none(value: str) -> bool:
    if value is None:
        return True
    return False

def is_int(component: Optional[InputComponent | OutputComponent], value: str) -> bool:
    # if component and not component.array and component.datatype.classname == 'Int':
    #     return True
    if expressions.is_int(str(value)):
        return True
    return False

def is_float(component: Optional[InputComponent | OutputComponent], value: str) -> bool:
    # if component and not component.array and component.datatype.classname == 'Float':
    #     return True
    if expressions.is_float(str(value)):
        return True
    return False
    
    

@dataclass
class StaticInputValue(InputValue):
    str_value: str
    is_default: bool

    def __post_init__(self):
        self.scatter: bool = False
    
    @property
    def value_type(self) -> str:
        """
        only StaticValueLinkingStrategy and DefaultValueLinkingStrategy 
        call select_input_value_type(). don't need to worry about CONNECTION and RUNTIME_VALUE
        """
        if is_bool(self.str_value):
            return 'boolean'
        elif is_int(self.component, self.str_value):
            return 'int'
        elif is_float(self.component, self.str_value):
            return 'float'
        elif is_none(self.str_value):
            return 'none'
        elif expressions.is_var(self.str_value) or expressions.has_var(self.str_value):
            return 'env_var'
        else:
            return 'string'
    
    @property
    def raw_value(self) -> Any:
        """
        get the typed value for this InputValue.
        literal_eval will take care of most cases.
        that said, sometimes it doesn't work for values which could be 
        interpreted as strings OR numeric types.

        eg self.string_value == '7', it will return as str('7'). 
        we may know this should actually be an int given other data. 
        """
        # None, True / False, string will be taken care of here.
        # try: except needed because literal_eval() will fail on specific strings
        val = self.str_value

        # empty string
        if val == '':       
            return None
        
        # list, tuple, dict
        for open_b, close_b in [('(', ')'), ('[', ']'), ('{', '}')]:
            if val[0] == open_b and val[-1] == close_b:
                val = ast.literal_eval(val)
        
        if val in ['True', 'False', 'None']:
            val = ast.literal_eval(val)
        
        # cast str to int or float if appropriate.
        if isinstance(val, str):
            if is_float(self.component, val):
                val = float(val)
            elif is_int(self.component, val):
                val = int(val)

        return val

    @property
    def wrapped_value(self) -> str:
        if self._should_wrap_value():
            return f'"{self.raw_value}"'
        else:
            return f'{self.raw_value}'
        
    def _should_wrap_value(self) -> bool:
        if self.value_type == 'string':
            return True
        if self.value_type == 'env_var':
            return True
        return False



@dataclass
class ConnectionInputValue(InputValue):
    step_uuid: str
    out_uuid: str
    
    @property
    def wrapped_value(self) -> str:
        step_tag = tags.get(self.step_uuid)
        out_tag = tags.get(self.out_uuid)
        return f'w.{step_tag}.{out_tag}'
    

@dataclass
class WorkflowInputInputValue(InputValue):
    input_uuid: str
    is_runtime: bool

    @property
    def wrapped_value(self) -> str:
        wflow_inp_tag = tags.get(self.input_uuid)
        return f'w.{wflow_inp_tag}'

