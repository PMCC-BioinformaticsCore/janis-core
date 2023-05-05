


from typing import Any
from copy import deepcopy
from janis_core.ingestion.galaxy.gx.gxtool.param import BoolParam
from janis_core.ingestion.galaxy.gx.gxtool import XMLToolDefinition


def resolve_values(gxstep: dict[str, Any], xmltool: XMLToolDefinition) -> dict[str, Any]:
    flattener = ToolStateValueResolver(xmltool)
    return flattener.flatten(gxstep['tool_state'])


class ToolStateValueResolver:
    def __init__(self, xmltool: XMLToolDefinition):
        self.xmltool = xmltool

    def flatten(self, the_dict: dict[str, Any]) -> dict[str, Any]:
        curr_path: list[str] = []
        input_dict = self.explore_node(the_dict, curr_path)
        return input_dict

    def explore_node(self, the_dict: dict[str, Any], path: list[str]) -> dict[str, Any]:
        keys_to_delete: list[str] = []

        for key, value in the_dict.items():
            if value is None:
                varname = self.get_path_as_str(key, path)
                varname = f'${varname}'
                the_dict[key] = varname

            elif value == {"__class__": "RuntimeValue"}:
                varname = self.get_path_as_str(key, path)
                varname = f'${varname}'
                the_dict[key] = varname
                # keys_to_delete.append(key)
                
            elif value == {"__class__": "ConnectedValue"}:
                varname = self.get_path_as_str(key, path)
                varname = f'${varname}'
                the_dict[key] = varname
                # keys_to_delete.append(key)
            
            # if value == {"__class__": "RuntimeValue"}:
            #     the_dict[key] = '__RuntimeValue__'
                
            # elif value == {"__class__": "ConnectedValue"}:
            #     the_dict[key] = '__ConnectedValue__'

            elif key == '__current_case__':
                keys_to_delete.append(key)

            elif isinstance(value, dict):
                curr_path = deepcopy(path)
                curr_path.append(key)
                the_dict[key] = self.explore_node(value, curr_path)  # type: ignore
            
            else:
                # resolving bool flag values
                gxvarname = self.get_path_as_str(key, path)
                param = self.xmltool.inputs.get(gxvarname)
                if param and isinstance(param, BoolParam):
                    if value == 'false':
                        the_dict[key] = param.falsevalue
                    else:
                        the_dict[key] = param.truevalue
                    continue
        
        for key in keys_to_delete:
            del the_dict[key]
        
        return the_dict

    def get_path_as_str(self, name: str, path_copy: list[str]) -> str:
        if len(path_copy) > 0:
            full_name = f'{".".join(path_copy)}.{name}'
        else:
            full_name = name
        return full_name
    