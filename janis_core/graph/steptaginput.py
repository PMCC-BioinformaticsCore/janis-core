from typing import Optional, Dict, Any

from janis_core.graph.node import Node, NodeTypes
from janis_core.tool.tool import ToolOutput, ToolInput
from janis_core.types.common_data_types import Array
from janis_core.utils import first_value
from janis_core.utils.logger import Logger


def full_lbl(node: Node, tag: Optional[str]) -> str:
    if tag is None:
        return node.id()
    return f"{node.id()}/{tag}"


def full_dot(node: Node, tag: Optional[str]) -> str:
    if tag is None:
        return node.id()
    return f"{node.id()}.{tag}"


class Edge:
    def __init__(
        self,
        start: Node,
        stag: Optional[str],
        finish: Node,
        ftag: Optional[str],
        should_scatter,
    ):
        Logger.log(
            f"Creating edge: ({NodeTypes.to_str(start.node_type)}) '{start.id()}.{stag}' → "
            f"({NodeTypes.to_str(finish.node_type)}) '{finish.id()}.{ftag}'"
        )

        self.start: Node = start
        self.stag: Optional[str] = stag
        self.finish: Node = finish
        self.ftag: Optional[str] = ftag
        self.compatible_types: Optional[bool] = None
        self.scatter = should_scatter

        self.validate_tags()
        self.check_types()

    def source_slashed(self):
        return full_lbl(self.start, self.stag)

    def source_dotted(self):
        return full_dot(self.start, self.stag)

    def validate_tags(self):
        if (
            self.start.node_type == NodeTypes.STEP
            and self.stag not in self.start.outputs()
        ):
            raise Exception(
                f"Could not find the tag '{self.stag}' in the inputs of '{self.start.id()}'"
            )
        if (
            self.finish.node_type == NodeTypes.STEP
            and self.ftag not in self.finish.inputs()
        ):
            raise Exception(
                f"Could not find the tag '{self.ftag}' in the outputs of '{self.finish.id()}': {list(self.finish.inputs().keys())}"
            )

    def check_types(self):
        from janis_core.workflow.workflow import InputNode, StepNode

        stoolin: ToolOutput = self.start.outputs()[
            self.stag
        ] if self.stag is not None else first_value(self.start.outputs())
        ftoolin: ToolInput = self.finish.inputs()[
            self.ftag
        ] if self.ftag is not None else first_value(self.finish.inputs())

        stype = stoolin.output_type
        ftype = ftoolin.input_type

        start_is_scattered = (
            isinstance(self.start, StepNode) and self.start.scatter is not None
        )

        if start_is_scattered:
            Logger.log(
                f"This edge merges the inputs from '{full_dot(self.start, self.stag)}' for "
                f"'{full_dot(self.finish, self.ftag)}'"
            )
            stype = Array(stype)

        if self.scatter:
            if not isinstance(stype, Array):
                raise Exception(
                    f"Scatter was required for '{self.start.id()}.{self.stag} → '{self.finish.id()}.{self.ftag}' but "
                    f"the input type was {type(stype).__name__} and not an array"
                )
            stype = stype.subtype()

        source_has_default = (
            isinstance(self.start, InputNode) and self.start.default is not None
        )

        # Scatters are handled automatically by the StepTagInput Array unwrapping
        # Merges are handled automatically by the `start_is_scattered` Array wrap

        self.compatible_types = ftype.can_receive_from(stype, source_has_default)
        if not self.compatible_types:
            if isinstance(ftype, Array) and ftype.subtype().can_receive_from(stype):
                self.compatible_types = True

        if not self.compatible_types:

            s = full_dot(self.start, self.stag)
            f = full_dot(self.finish, self.ftag)
            message = (
                f"Mismatch of types when joining '{s}' to '{f}': "
                f"{stoolin.output_type.id()} -/→ {ftoolin.input_type.id()}"
            )
            if isinstance(stype, Array) and ftype.can_receive_from(stype.subtype()):
                message += " (did you forget to SCATTER?)"
            Logger.critical(message)


class StepTagInput:
    """
    This class represents the connections that a single input on a step has. Hence, a step
    will have one StepTagInput for each potential input of the tool.
    """

    def __init__(self, finish: Node, finish_tag: str):

        self.finish: Node = finish
        self.ftag: Optional[str] = finish_tag

        self.default = None
        self.multiple_inputs = False

        self.source_map: Dict[str, Edge] = {}

    def add_source(self, start: Node, stag: Optional[str], should_scatter) -> Edge:
        """
        Add a connection
        :param start:
        :param stag:
        :param should_scatter:
        :return:
        """

        from janis_core.workflow.workflow import StepNode

        stype = (
            start.outputs()[stag] if stag is not None else first_value(start.outputs())
        ).output_type
        ftype = (
            self.finish.inputs()[self.ftag]
            if self.ftag is not None
            else first_value(self.finish.inputs())
        ).input_type

        start_is_scattered = isinstance(start, StepNode) and start.scatter is not None

        if start_is_scattered:
            Logger.log(
                f"This edge merges the inputs from '{full_dot(start, stag)}' for "
                f"'{full_dot(self.finish, self.ftag)}'"
            )
            stype = Array(stype)

        if should_scatter:
            if not isinstance(stype, Array):
                raise Exception(
                    f"Scatter was required for '{start.id()}.{stag} → '{self.finish.id()}.{self.ftag}' but "
                    f"the input type was {type(stype).__name__} and not an array"
                )
            stype = stype.subtype()

        if len(self.source_map) == 1 and start.id() not in self.source_map:
            self.multiple_inputs = True

            if not isinstance(ftype, Array):
                raise Exception(
                    f"Adding multiple inputs to '{self.finish.id()}' and '{ftype.id()}' is not an array"
                )

        if not isinstance(stype, Array) and isinstance(ftype, Array):
            # https://www.commonwl.org/user_guide/misc/#connect-a-solo-value-to-an-input-that-expects-an-array-of-that-type
            self.multiple_inputs = True

        e = Edge(start, stag, self.finish, self.ftag, should_scatter=should_scatter)
        self.source_map[start.id()] = e
        return e

    def set_default(self, default: Any):
        Logger.log(
            f"Setting the default of '{self.finish.id()}.{self.ftag}' to be '{str(default)}'"
        )
        self.default = default

    def source(self):
        n = len(self.source_map)
        if n == 0:
            return None
        elif n == 1:
            return first_value(self.source_map)
        else:
            return list(self.source_map.values())

    def dotted_source(self):
        n = len(self.source_map)

        if n == 0:
            return None
        elif n == 1:
            return first_value(self.source_map).source_dotted()
        else:
            return [e.source_dotted() for e in self.source_map.values()]

    def slashed_source(self):
        n = len(self.source_map)
        if n == 0:
            return None
        elif n == 1:
            return first_value(self.source_map).source_slashed()
        else:
            return [e.source_slashed() for e in self.source_map.values()]
