

from typing import Optional
from dataclasses import dataclass, field

from janis_core import (
    ToolInput,
    TInput,
)
from ..scope import Scope


@dataclass
class ProcessDSCategoryRegister:
    
    data_structure: dict[str, dict[str, set[str]]] = field(default_factory=dict)
    
    """
    for each CommandTool / PythonTool, stores data on which ToolInputs are:
        - process inputs
        - param inputs
        - internal inputs

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
class ProcessDSVariableRegister:
    
    data_structure: dict[str, dict[str, Optional[str | list[str]]]] = field(default_factory=dict)
    
    """
    for each CommandTool / PythonTool, stores a map of tinput_id: variable_name.
    the variable_name is how the particular tinput_id will be referenced inside the process. 
    
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

    

pds_register = ProcessDSCategoryRegister()
pvn_register = ProcessDSVariableRegister()

# categories
def update_categories(scope: Scope, process_inputs: set[str], param_inputs: set[str], internal_inputs: set[str]):
    pds_register.add(scope, 'process', process_inputs)
    pds_register.add(scope, 'param', param_inputs)
    pds_register.add(scope, 'internal', internal_inputs)

def process_inputs(scope: Scope) -> set[str]:
    return pds_register.get(scope, 'process')

def param_inputs(scope: Scope) -> set[str]:
    return pds_register.get(scope, 'param')

def internal_inputs(scope: Scope) -> set[str]:
    return pds_register.get(scope, 'internal')

# variable references
def update_variables(scope: Scope, tinput_id: str, name: Optional[str | list[str]]):
    pvn_register.add(scope, tinput_id, name)

def get_variable(scope: Scope, inp: ToolInput | TInput) -> Optional[str | list[str]]:
    return pvn_register.get(scope, inp.id())