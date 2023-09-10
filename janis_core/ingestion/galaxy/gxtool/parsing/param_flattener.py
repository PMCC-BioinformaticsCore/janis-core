


from copy import deepcopy
from typing import Any

from galaxy.tools.parameters.basic import ToolParameter, HiddenToolParameter
from galaxy.tools.parameters.grouping import Conditional, ConditionalWhen, Section, Repeat

XmlNode = ToolParameter | Conditional | ConditionalWhen | Section | Repeat


class XMLParamFlattener:
    def __init__(self, inputs: dict[str, Any]):
        self.inputs = inputs
        self.flat_params: list[ToolParameter] = [] 

    def flatten(self) -> list[ToolParameter]:
        for node in self.inputs.values():
            self.explore_node(node, [])
        return self.flat_params

    def explore_node(self, node: XmlNode, heirarchy: list[str]) -> None:
        heirarchy = deepcopy(heirarchy)
        # TODO move these to check attributes to remove galaxy class imports
        match node:
            case HiddenToolParameter():
                pass
            case ToolParameter(): # leaf node
                self.flatten_param(node, heirarchy)
            case Conditional(): 
                heirarchy.append(node.name)
                self.flatten_param(node.test_param, heirarchy)
                for child in node.cases:
                    self.explore_node(child, heirarchy)
            case ConditionalWhen():
                for child in node.inputs.values():
                    self.explore_node(child, heirarchy)
            case Section():
                heirarchy.append(node.name)
                for child in node.inputs.values():
                    self.explore_node(child, heirarchy)
            case Repeat():
                pass
                # logging.has_repeat()
            case _:
                raise NotImplementedError()
    
    def flatten_param(self, param: ToolParameter, heirarchy: list[str]) -> None:
        prefix = ".".join(heirarchy) + '.' if len(heirarchy) > 0 else ''
        flat_name = f'{prefix}{param.name}'
        param.flat_name = flat_name
        self.flat_params.append(param)

