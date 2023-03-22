


from dataclasses import dataclass, field
from ...scope import Scope


@dataclass
class ProcessDataSourceRegister:
    
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
    


pds_register = ProcessDataSourceRegister()

def update(scope: Scope, process_inputs: set[str], param_inputs: set[str], internal_inputs: set[str]):
    pds_register.add(scope, 'process', process_inputs)
    pds_register.add(scope, 'param', param_inputs)
    pds_register.add(scope, 'internal', internal_inputs)

def process_inputs(scope: Scope) -> set[str]:
    return pds_register.get(scope, 'process')

def param_inputs(scope: Scope) -> set[str]:
    return pds_register.get(scope, 'param')

def internal_inputs(scope: Scope) -> set[str]:
    return pds_register.get(scope, 'internal')


