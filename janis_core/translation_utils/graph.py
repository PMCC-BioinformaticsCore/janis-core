

from typing import Any
from janis_core import TInput
from janis_core.types import DataType
from janis_core.workflow.workflow import InputNode
from janis_core.graph.steptaginput import StepTagInput
from janis_core.operators.selectors import InputNodeSelector, StepOutputSelector

from .datatypes import get_dtt, DTypeType


def resolve_node(node: Any) -> Any:
    if isinstance(node, StepTagInput):
        return resolve_node(node.source_map[0].source)
    elif isinstance(node, InputNodeSelector):
        return node.input_node
    elif isinstance(node, StepOutputSelector):
        return node.node
    else:
        return node
    
def looks_like_placeholder_node(node: InputNode, step_id: str, tinput_id: str, tinput_dtype: DataType) -> bool:
    # single reference with dummy format name implies placeholder
    # keep file types
    dtt = get_dtt(node.datatype)
    if dtt in [
        DTypeType.SECONDARY_ARRAY,
        DTypeType.SECONDARY,
        DTypeType.FILE_PAIR_ARRAY,
        DTypeType.FILE_PAIR,
        DTypeType.FILE_ARRAY,
        DTypeType.FILE,
    ]:
        return False
    
    # node must have default
    if node.default is None:
        return False
    
    # name must look like placeholder
    if node.id() != f'{step_id}_{tinput_id}':
        return False
    
    # datatypes must match
    if node.datatype.__class__.__name__ != tinput_dtype.__class__.__name__:
        return True
    
    return True