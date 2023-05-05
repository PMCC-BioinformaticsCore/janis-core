


from typing import Any
from copy import deepcopy


def get_flattened_tool_state_old(gxstep: dict[str, Any]) -> dict[str, Any]:
    flattener = ToolStateFlattener()
    return flattener.flatten(gxstep)


class ToolStateFlattener:
    def __init__(self):
        self.flattened_tool_state: dict[str, Any] = {}

    def flatten(self, step: dict[str, Any]) -> dict[str, Any]:
        for name, value in step['tool_state'].items():
            self.explore_node(name, value, [])
        return self.flattened_tool_state

    def explore_node(self, name: str, value: Any, path: list[str]) -> Any:
        path_copy = deepcopy(path)
        if name == '__current_case__':
            pass
        elif value == {"__class__": "RuntimeValue"}:
            self.add_to_flattened_tool_state(name, '__RuntimeValue__', path_copy)
        elif value == {"__class__": "ConnectedValue"}:
            pass
        elif isinstance(value, dict):
            path_copy.append(name)
            for key, val in value.items():
                self.explore_node(key, val, path_copy)
        else:
            self.add_to_flattened_tool_state(name, value, path_copy)
    
    def add_to_flattened_tool_state(self, name: str, value: Any, path_copy: list[str]) -> None:
        if len(path_copy) > 0:
            full_name = f'{".".join(path_copy)}.{name}'
        else:
            full_name = name
        self.flattened_tool_state[full_name] = value


