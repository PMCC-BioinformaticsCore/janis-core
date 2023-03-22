

from typing import Optional
from dataclasses import dataclass, field

from janis_core import (
    ToolInput,
    TInput,
)

from ..scope import Scope



@dataclass
class ProcessVariableNameRegister:
    
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
    

pvn_register = ProcessVariableNameRegister()

def update(scope: Scope, tinput_id: str, name: Optional[str | list[str]]):
    pvn_register.add(scope, tinput_id, name)

def get(scope: Scope, inp: ToolInput | TInput) -> Optional[str | list[str]]:
    return pvn_register.get(scope, inp.id())

