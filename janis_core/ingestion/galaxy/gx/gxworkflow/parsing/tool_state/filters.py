




from typing import Any
from abc import ABC, abstractmethod
from copy import deepcopy
from janis_core.ingestion.galaxy.gx.gxtool.param import BoolParam
from janis_core.ingestion.galaxy.gx.gxtool import XMLToolDefinition


### helper functions ###

def get_path_as_str(name: str, path_copy: list[str]) -> str:
    if len(path_copy) > 0:
        full_name = f'{".".join(path_copy)}.{name}'
    else:
        full_name = name
    return full_name


### filters ###

class Filter(ABC):
    def __init__(self, the_dict: dict[str, Any], path: list[str], xmltool: XMLToolDefinition) -> None:
        self.the_dict = the_dict
        self.path = path
        self.xmltool = xmltool
    
    @abstractmethod
    def apply(self) -> dict[str, Any]:
        ...


class ReplaceNullWithVarname(Filter):
    
    def apply(self) -> dict[str, Any]:
        for key, value in self.the_dict.items():
            if value is None:
                varname = get_path_as_str(key, self.path)
                varname = f'${varname}'
                self.the_dict[key] = varname
        return self.the_dict


class ReplaceConnectedWithVarname(Filter):
    
    def apply(self) -> dict[str, Any]:
        for key, value in self.the_dict.items():
            if value == {"__class__": "ConnectedValue"}:
                varname = get_path_as_str(key, self.path)
                varname = f'${varname}'
                self.the_dict[key] = varname
        return self.the_dict


class DeNestClass(Filter):
    
    def apply(self) -> dict[str, Any]:
        for key, value in self.the_dict.items():
            if value == {"__class__": "ConnectedValue"} or value == {"__class__": "RuntimeValue"}:
                self.the_dict[key] = value["__class__"]
        return self.the_dict


class ReplaceRuntimeWithVarname(Filter):
    
    def apply(self) -> dict[str, Any]:
        for key, value in self.the_dict.items():
            if value == {"__class__": "RuntimeValue"}:
                varname = get_path_as_str(key, self.path)
                varname = f'${varname}'
                self.the_dict[key] = varname
        return self.the_dict


class IgnoreCurrentCase(Filter):
    
    def apply(self) -> dict[str, Any]:
        keys_to_delete: list[str] = []
        for key in self.the_dict.keys():
            if key == '__current_case__':
                keys_to_delete.append(key)
        
        for key in keys_to_delete:
            del self.the_dict[key]
        return self.the_dict


class ReplaceBoolWithValue(Filter):
    
    def apply(self) -> dict[str, Any]:
        for key, value in self.the_dict.items():
            gxvarname = get_path_as_str(key, self.path)
            param = self.xmltool.inputs.get(gxvarname)
            if param and isinstance(param, BoolParam):
                if value == 'false':
                    self.the_dict[key] = param.falsevalue
                else:
                    self.the_dict[key] = param.truevalue
        return self.the_dict


class Flatten(ABC):
    
    def __init__(self, the_dict: dict[str, Any], path: list[str], xmltool: XMLToolDefinition) -> None:
        self.the_dict = the_dict
        self.path = path
        self.xmltool = xmltool
        self.flattened_tool_state: dict[str, Any] = {}
    
    def apply(self) -> dict[str, Any]:
        curr_path: list[str] = []
        self.explore_node(self.the_dict, curr_path)
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

