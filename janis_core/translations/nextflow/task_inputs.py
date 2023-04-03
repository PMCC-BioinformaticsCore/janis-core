

from dataclasses import dataclass, field

@dataclass
class MinimalTaskInputRegister:
    """for each task to call, stores minimal tool inputs needed as task inputs."""
    data_structure: dict[str, set[str]] = field(default_factory=dict)
    
    def to_string(self) -> str:
        out: str = ''
        for tool_id, minimal_task_input_ids in self.data_structure.items():
            out += f'\n{tool_id} ---\n'
            for tinput_id in minimal_task_input_ids:
                out += f'{tinput_id}\n'
        return out

mti_register = MinimalTaskInputRegister()

def add(tool_id: str, minimal_task_input_ids: set[str]) -> None:
    mti_register.data_structure[tool_id] = minimal_task_input_ids

def get(tool_id: str) -> set[str]:
    return mti_register.data_structure[tool_id]

def clear() -> None:
    mti_register.data_structure = {}    