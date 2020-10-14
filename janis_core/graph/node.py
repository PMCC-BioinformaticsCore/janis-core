"""
    node.py

    Provides base class that different nodes must override, this translates closest to a Step
"""
from abc import ABC, abstractmethod
from enum import Enum
from typing import Dict, List, Tuple, Any

from janis_core.tool.tool import TInput, TOutput

NodeLabel = str


class NodeType(Enum):
    INPUT = 1
    OUTPUT = 2
    STEP = 3

    @staticmethod
    def to_str(node_type) -> str:

        v = NodeType(node_type) if not isinstance(node_type, NodeType) else node_type

        if v == NodeType.INPUT:
            return "Input"
        elif v == NodeType.OUTPUT:
            return "Output"
        elif v == NodeType.STEP:
            return "Task"
        raise Exception(f"Unhandled task type: '{v.value}'")

    def __str__(self):
        return NodeType.to_str(self)

    @staticmethod
    def to_col(node_type) -> str:
        if node_type == NodeType.INPUT:
            return "red"
        if node_type == NodeType.OUTPUT:
            return "lightblue"
        if node_type == NodeType.STEP:
            return "blue"
        raise Exception(f"Unhandled task type: '{node_type}'")


class Node(ABC):

    _N_counter: int = 1
    _N_nodeId_map: Dict[int, Any] = {}

    def __init__(self, wf, node_type: NodeType, identifier: NodeLabel, depth=0):

        self.wf = wf
        self.node_type: NodeType = node_type
        self.identifier: NodeLabel = identifier
        self.depth = depth

        # actually a StepInput
        self.sources: Dict[str, any] = {}

        # Update unique counter for hash
        self._nodeId = Node._N_counter
        Node._N_counter += 1

        # Map the node, so we can look it up later
        self._N_nodeId_map[self._nodeId] = self

    def id(self) -> str:
        return self.identifier

    def __hash__(self):
        return self._nodeId

    def __repr__(self):
        return f"{self.node_type}: {self.id()}"

    def __str__(self):
        return f"{NodeType.to_str(self.node_type)}: {self.identifier}"

    def set_depth(self, depth: int):
        self.depth = max(self.depth, depth)

    def __eq__(self, other):
        if isinstance(other, Node):
            return self._nodeId == other._nodeId
        return False

    @abstractmethod
    def inputs(self) -> Dict[str, TInput]:
        raise Exception(
            f"Subclass {type(self)} must implement inputs, return dict: key: ToolInput"
        )

    @abstractmethod
    def outputs(self) -> Dict[str, TOutput]:
        raise Exception(
            f"Subclass {type(self)} must implement outputs, return dict: key: ToolOutput"
        )


NodeAndTag = Tuple[str, Node]


def layout_nodes(nodes: List[Node], n_inputs: int = 0) -> Dict[Node, Tuple[int, int]]:
    """
    Stack on depth away from root, and scale smaller columns
    :param n_inputs:
    :param nodes: list of nodes in Graph to position
    :return: Dict[Node, (x,y)]
    """
    pos = {}
    # int is the depth at the specified index
    cur_depth: List[int] = []
    depth_node: Dict[int, List[Node]] = {}

    def get_idx_for_depth(d: int):
        d_idx = d
        m = len(cur_depth)
        if d_idx >= m:
            cur_depth.extend([0] * (d - m + 1))
        return cur_depth[d]

    def push_to_dict(key, val):
        if key in depth_node:
            depth_node[key].append(val)
        else:
            depth_node[key] = [val]

    for node in nodes:
        depth = node.depth
        push_to_dict(depth, node)
        d = get_idx_for_depth(depth)
        cur_depth[depth] += 1
        pos[node] = (depth, d)

    # Now normalise each depth
    max_in_col = float(max(cur_depth + [n_inputs]))

    for (idx, d) in enumerate(cur_depth):
        if idx not in depth_node:
            continue
        nodes = depth_node[idx]
        scale = max_in_col / len(nodes)

        # If max and this col differ in even / odd, add half a unit
        half_bias = 0 if (len(nodes) % 2) == (max_in_col % 2) else 0.5

        for node in nodes:
            (x, y) = pos[node]
            pos[node] = (x, (y + half_bias) * scale)

    return pos


def layout_nodes2(nodes: List[Node]) -> Dict[Node, Tuple[int, int]]:
    # Aim for something like Rabix
    inputs = [n for n in nodes if n.node_type == NodeType.INPUT]
    others = [n for n in nodes if n.node_type != NodeType.INPUT]

    pos = layout_nodes(others, len(inputs))
    s = 0
    for n in inputs:
        pos[n] = (0, s)
        s += 1
    return pos


def layout_nodes3(nodes: List[Node]) -> Dict[Node, Tuple[int, int]]:
    # Other ideas to scale might be, start with steps and
    # stack in a similar way to depth around center, then do the same
    return {}
