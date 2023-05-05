



from typing import Any
from copy import deepcopy


def flatten(inputs_dict: dict[str, Any]) -> dict[str, Any]:
    flattener = ToolStateFlattener()
    return flattener.flatten(inputs_dict)


class ToolStateFlattener:
    def __init__(self):
        self.flattened_tool_state: dict[str, Any] = {}

    def flatten(self, the_dict: dict[str, Any]) -> dict[str, Any]:
        curr_path: list[str] = []
        self.explore_node(the_dict, curr_path)
        return self.flattened_tool_state

    def explore_node(self, the_dict: dict[str, Any], path: list[str]) -> None:
        for key, value in the_dict.items():
            if isinstance(value, dict):
                curr_path = deepcopy(path)
                curr_path.append(key)
                self.explore_node(value, curr_path)  # type: ignore
            else:
                self.add_to_flattened_tool_state(key, value, path)
    
    def add_to_flattened_tool_state(self, name: str, value: Any, path_copy: list[str]) -> None:
        if len(path_copy) > 0:
            full_name = f'{".".join(path_copy)}.{name}'
        else:
            full_name = name
        self.flattened_tool_state[full_name] = value
