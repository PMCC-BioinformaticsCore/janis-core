from typing import Optional

from janis_core.types import get_instantiated_type
from janis_core.operators import Selector
from janis_core.operators import Operator
from janis_core.graph.node import Node, NodeType
from janis_core.tool.tool import TInput
from janis_core.utils import first_value
from janis_core.utils.logger import Logger
from janis_core import settings
from janis_core.messages import log_message
from janis_core.messages import ErrorCategory
from uuid import uuid4


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
        self, source: Selector, finish: Node, ftag: Optional[str], should_scatter: Optional[bool]=None
    ):
        Logger.log(
            f"Creating edge: ({source} → "
            f"({NodeType.to_str(finish.node_type)}) '{finish.id()}.{ftag}'"
        )
        self.uuid: str = str(uuid4())
        self.source = source
        self.finish: Node = finish
        self.ftag: Optional[str] = ftag
        self.compatible_types: Optional[bool] = None
        self.should_scatter = should_scatter
        self.validate_tags()
        self.check_types()

    def validate_tags(self):
        if self.finish.node_type == NodeType.STEP:
            if self.ftag not in self.finish.inputs():
                if settings.graph.ALLOW_UNKNOWN_SOURCE:
                    msg = "Could not connect this input to its data source"
                    log_message(self.uuid, msg, ErrorCategory.PLUMBING)
                else:
                    raise Exception(
                        f"Could not find the tag '{self.ftag}' in the outputs of '{self.finish.id()}': {list(self.finish.inputs().keys())}"
                    )

    def check_types(self):
        # stoolin: TOutput = self.start.outputs()[
        #     self.stag
        # ] if self.stag is not None else first_value(self.start.outputs())

        ftoolin: TInput = (
            self.finish.inputs()[self.ftag]
            if self.ftag is not None
            else first_value(self.finish.inputs())
        )

        stype = get_instantiated_type(self.source.returntype())
        ftype = get_instantiated_type(ftoolin.intype)

        if self.should_scatter:
            if not stype.is_array():
                if settings.graph.ALLOW_NON_ARRAY_SCATTER_INPUT:
                    msg = f"This task is supposed to run in parallel across this input ({ftoolin.id()}), but the data source is not an array."
                    log_message(self.uuid, msg, ErrorCategory.PLUMBING)
                else:
                    raise Exception(
                        f"Scatter was required for '{operator} → '{self.finish.id()}.{self.ftag}' but "
                        f"the input type was {type(stype).__name__} and not an array"
                    )
            else:
                stype = stype.subtype()

        # Scatters are handled automatically by the StepTagInput Array unwrapping
        # Merges are handled automatically by the `start_is_scattered` Array wrap

        self.compatible_types = ftype.can_receive_from(stype, False)
        if not self.compatible_types:
            if ftype.is_array() and ftype.subtype().can_receive_from(stype):
                self.compatible_types = True

        if not self.compatible_types:
            if settings.graph.ALLOW_INCOMPATIBLE_TYPES:
                msg = f"The data source for this input is a {stype.id()}, but the input is a {ftoolin.intype.id()}"
                log_message(self.uuid, msg, ErrorCategory.DATATYPES)
            else:
                s = str(self.source)
                f = full_dot(self.finish, self.ftag)
                message = (
                    f"Mismatch of types when joining '{s}' to '{f}': "
                    f"{stype.id()} -/→ {ftoolin.intype.id()}"
                )
                if stype.is_array() and ftype.can_receive_from(stype.subtype()):
                    message += " (did you forget to SCATTER?)"
                Logger.critical(message)


class StepTagInput:
    """
    This class represents the connections that a single input on a step has. 
    A step will have one StepTagInput for each potential input of the tool.
    """

    def __init__(self, finish: Node, finish_tag: str):
        self.uuid: str = str(uuid4())
        self.finish: Node = finish
        self.ftag: Optional[str] = finish_tag
        self.multiple_inputs = False
        self.source_map: list[Edge] = []
        self.operator: Optional[Operator] = None

    def add_source(self, operator: Selector, should_scatter: Optional[bool]=None) -> Edge:
        """
        Add a connection
        :param start:
        :param stag:
        :param should_scatter:
        :return:
        """
        stype = get_instantiated_type(operator.returntype())

        if self.ftag:
            tinput = self.finish.inputs()[self.ftag]
        else:
            tinput = first_value(self.finish.inputs())        
        ftype = tinput.intype

        # from janis_core.workflow.workflow import StepNode
        # start_is_scattered = isinstance(start, StepNode) and start.scatter is not None
        #
        # if start_is_scattered:
        #     Logger.log(
        #         f"This edge merges the inputs from '{full_dot(start, stag)}' for "
        #         f"'{full_dot(self.finish, self.ftag)}'"
        #     )
        #     stype = Array(stype)

        if should_scatter:
            if not stype.is_array():
                if settings.graph.ALLOW_NON_ARRAY_SCATTER_INPUT:
                    msg = f"This task is supposed to run in parallel across this input ({tinput.id()}), but the data source is not an array."
                    log_message(self.uuid, msg, ErrorCategory.PLUMBING)
                else:
                    raise Exception(
                        f"Scatter was required for '{operator} → '{self.finish.id()}.{self.ftag}' but "
                        f"the input type was {type(stype).__name__} and not an array"
                    )
            else:
                stype = get_instantiated_type(stype.subtype())

        if len(self.source_map) == 1:  # and start.id() not in self.source_map:
            self.multiple_inputs = True
            if not ftype.is_array():
                if settings.graph.ALLOW_INCORRECT_NUMBER_OF_SOURCES:
                    msg = "This input has multiple data sources, but should only have one (it is not an array)."
                    log_message(self.uuid, msg, ErrorCategory.PLUMBING)
                else:
                    raise Exception(
                        f"Adding multiple inputs to '{self.finish.id()}' "
                        f"and '{ftype.id()}' is not an array"
                    )

        if not stype.is_array() and ftype.is_array():
            # https://www.commonwl.org/user_guide/misc/#connect-a-solo-value-to-an-input-that-expects-an-array-of-that-type
            self.multiple_inputs = True

        e = Edge(operator, self.finish, self.ftag, should_scatter=should_scatter)
        # todo: deal with source_map
        self.source_map.append(e)
        return e

    def source(self):
        n = len(self.source_map)
        if n == 0:
            return None
        elif n == 1:
            return self.source_map[0]
        else:
            return list(self.source_map)
