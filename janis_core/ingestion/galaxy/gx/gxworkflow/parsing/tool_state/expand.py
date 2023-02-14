



import json
from typing import Any


def expand_tool_state(step: dict[str, Any]) -> dict[str, Any]:
    expander = ToolStateExpander(step)
    return expander.load()


class ToolStateExpander:
    def __init__(self, step: dict[str, Any]):
        self.step = step
        self.tool_state: dict[str, Any] = {}
        #self.tool_state = json.loads(step['tool_state'])

    def load(self) -> dict[str, Any]:
        tool_state = self.explore_node(self.step['tool_state'])
        return tool_state
    
    def explore_node(self, node: Any) -> Any:
        if self.is_quoted(node) or self.is_brackets(node):
            node = json.loads(node)
            if isinstance(node, dict):
                for key, val in node.items():
                    node[key] = self.explore_node(val) # type: ignore
        return node

    def is_quoted(self, node: Any) -> bool:
        if isinstance(node, str):
            if node.startswith('"') and node.endswith('"'): 
                return True
            elif node.startswith("'") and node.endswith("'"): 
                return True
        return False
    
    def is_brackets(self, node: Any) -> bool:
        if isinstance(node, str):
            if node.startswith('{') and node.endswith('}'): 
                return True
        return False

    

