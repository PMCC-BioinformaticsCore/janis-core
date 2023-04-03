

from typing import Optional
from dataclasses import dataclass, field

from janis_core import (
    ToolInput,
    TInput,
)
from .scope import Scope


@dataclass
class TaskDSCategoryRegister:
    
    data_structure: dict[str, dict[str, set[str]]] = field(default_factory=dict)
    
    """
    for each Workflow / CommandTool / PythonTool, for each ToolInput, stores data on whether the ToolInput is fed value via:
        - a process input
        - a global param
        - not fed a value (internal input - static value)

    data_structure description:
    {
        Scope: {
            'process': set[str],
            'param': set[str],
            'internal': set[str]
        }
    }
    """
    
    def add(self, scope: Scope, subtype: str, tinput_ids: set[str]) -> None:
        """
        adds some data to our data_structure. 
        scope: tool scope
        subtype: one of 'process', 'param', 'internal'
        tinput_ids: set of ToolInput identifiers belonging to that subtype
        """
        label = scope.to_string()
        if label not in self.data_structure:
            self.data_structure[label] = {}
        self.data_structure[label][subtype] = tinput_ids

    def get(self, scope: Scope, subtype: str) -> set[str]:
        label = scope.to_string()
        return self.data_structure[label][subtype]

    def to_string(self) -> str:
        out: str = ''
        for scope_label, categories in self.data_structure.items():
            out += f'\n{scope_label}\n'
            for catname, tinput_ids in categories.items():
                out += f'{catname}: {" ".join(tinput_ids)}\n'
        return out
    


@dataclass
class TaskDSVariableRegister:
    
    data_structure: dict[str, dict[str, Optional[str | list[str]]]] = field(default_factory=dict)
    
    """ 
    for each Workflow / CommandTool / PythonTool, stores a map of tinput_id: data_source.
    the data_source is the initial value which drives the ToolInput inside a task. 
    
    data_structure description:
    {
        Scope: {
            'my_file': str,
            'my_string': None,
            'my_secondary': list[str],
            'my_secondary_arr': str,
        }
    }
    """
    
    def add(self, scope: Scope, tinput_id: str, name: Optional[str | list[str]]) -> None:
        """
        adds some data to our data_structure. 
        scope: tool scope
        
        """
        label = scope.to_string()
        if label not in self.data_structure:
            self.data_structure[label] = {}
        self.data_structure[label][tinput_id] = name

    def get(self, scope: Scope, tinput_id: str) -> Optional[str | list[str]]:
        label = scope.to_string()
        return self.data_structure[label][tinput_id]
    
    def to_string(self) -> str:
        out: str = ''
        for scope_label, name_map in self.data_structure.items():
            out += f'\n{scope_label}\n'
            for tinput_id, varname in name_map.items():
                out += f'{tinput_id}: {varname}\n'
        return out

    

tds_register = TaskDSCategoryRegister()
tvn_register = TaskDSVariableRegister()

# categories
def update_categories(scope: Scope, process_inputs: set[str], param_inputs: set[str], internal_inputs: set[str]):
    tds_register.add(scope, 'process', process_inputs)
    tds_register.add(scope, 'param', param_inputs)
    tds_register.add(scope, 'internal', internal_inputs)

def process_inputs(scope: Scope) -> set[str]:
    return tds_register.get(scope, 'process')

def param_inputs(scope: Scope) -> set[str]:
    return tds_register.get(scope, 'param')

def internal_inputs(scope: Scope) -> set[str]:
    return tds_register.get(scope, 'internal')

# variable references
def update_variables(scope: Scope, tinput_id: str, name: Optional[str | list[str]]):
    tvn_register.add(scope, tinput_id, name)

def get_variable(scope: Scope, inp: ToolInput | TInput) -> Optional[str | list[str]]:
    return tvn_register.get(scope, inp.id())

def clear() -> None:
    tds_register.data_structure = {}
    tvn_register.data_structure = {}