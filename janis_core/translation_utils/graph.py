

from typing import Any
from janis_core.graph.steptaginput import StepTagInput
from janis_core.operators.selectors import InputNodeSelector, StepOutputSelector


def resolve_node(node: Any) -> Any:
    if isinstance(node, StepTagInput):
        return resolve_node(node.source_map[0].source)
    elif isinstance(node, InputNodeSelector):
        return node.input_node
    elif isinstance(node, StepOutputSelector):
        return node.node
    else:
        return node