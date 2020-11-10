from abc import ABC, abstractmethod
from typing import Union

from janis_core.types import get_instantiated_type

from janis_core.graph.node import NodeType
from janis_core.types.common_data_types import Array, File, Directory, Int
from janis_core.utils import first_value
from janis_core.utils.logger import Logger


class Selector(ABC):
    @staticmethod
    def is_selector():
        return True

    @abstractmethod
    def returntype(self):
        pass

    def requires_contents(self):
        """
        A subclass should set this to TRUE
        :return:
        """
        return False

    def __neg__(self):
        self.negate()

    def negate(self):
        from janis_core.operators.logical import NotOperator

        return NotOperator(self)

    @abstractmethod
    def to_string_formatter(self):
        pass

    def assert_not_null(self):
        from janis_core.operators.logical import AssertNotNull

        return AssertNotNull(self)

    def is_null(self):
        return self.is_not_null().negate()

    def is_not_null(self):
        from janis_core.operators.logical import IsDefined

        return IsDefined(self)

    def and_(self, other):
        from janis_core.operators.logical import AndOperator

        return AndOperator(self, other)

    def __and__(self, other):
        return self.and_(other)

    def __rand__(self, other):
        from janis_core.operators.logical import AndOperator

        return AndOperator(other, self)

    def or_(self, other):
        from janis_core.operators.logical import OrOperator

        return OrOperator(self, other)

    def __or__(self, other):
        return self.or_(other)

    def __ror__(self, other):
        from janis_core.operators.logical import OrOperator

        return OrOperator(other, self)

    def add(self, other):
        from janis_core.operators.logical import AddOperator

        return AddOperator(self, other)

    def __add__(self, other):
        return self.add(other)

    def __radd__(self, other):
        from janis_core.operators.logical import AddOperator

        return AddOperator(other, self)

    def subtract(self, other):
        from janis_core.operators.logical import SubtractOperator

        return SubtractOperator(self, other)

    def __sub__(self, other):
        return self.subtract(other)

    def __rsub__(self, other):
        from janis_core.operators.logical import SubtractOperator

        return SubtractOperator(other, self)

    def multiply(self, other):
        from janis_core.operators.logical import MultiplyOperator

        return MultiplyOperator(self, other)

    def __mul__(self, other):
        return self.multiply(other)

    def __rmul__(self, other):
        from janis_core.operators.logical import MultiplyOperator

        return MultiplyOperator(other, self)

    def divide(self, other):

        from janis_core.operators.logical import DivideOperator

        return DivideOperator(self, other)

    def __truediv__(self, other):
        return self.divide(other)

    def __rtruediv__(self, other):
        from janis_core.operators.logical import DivideOperator

        return DivideOperator(other, self)

    # def __eq__(self, other):
    def equals(self, other):
        from janis_core.operators.logical import EqualityOperator

        return EqualityOperator(self, other)

    def not_equals(self, other):

        from janis_core.operators.logical import EqualityOperator

        return EqualityOperator(self, other)

    def __ne__(self, other):
        return self.not_equals(other)

    def greater_than(self, other):
        from janis_core.operators.logical import GtOperator

        return GtOperator(self, other)

    def __gt__(self, other):
        return self.greater_than(other)

    def greater_than_or_equals(self, other):
        from janis_core.operators.logical import GteOperator

        return GteOperator(self, other)

    def __ge__(self, other):
        return self.greater_than_or_equals(other)

    def less_than(self, other):
        from janis_core.operators.logical import LtOperator

        return LtOperator(self, other)

    def __lt__(self, other):
        return self.less_than(other)

    def less_than_or_equals(self, other):
        from janis_core.operators.logical import LteOperator

        return LteOperator(self, other)

    def length(self):
        from janis_core.operators.standard import LengthOperator

        return LengthOperator(self)

    def __len__(self):

        raise Exception(
            f"Calling the len function on a Janis selector, ie:'len({str(self)})' is not supported, please use '{str(self)}.length)'"
        )

    def as_str(self):
        from janis_core.operators.operator import AsStringOperator

        return AsStringOperator(self)

    def as_bool(self):
        from janis_core.operators.operator import AsBoolOperator

        return AsBoolOperator(self)

    def as_int(self):
        from janis_core.operators.operator import AsIntOperator

        return AsIntOperator(self)

    def as_float(self):
        from janis_core.operators.operator import AsFloatOperator

        return AsFloatOperator(self)

    def op_and(self, other):
        from janis_core.operators.logical import AndOperator

        return AndOperator(self, other)

    def op_or(self, other):
        from janis_core.operators.logical import OrOperator

        return OrOperator(self, other)

    def contents(self):
        from janis_core.operators.standard import ReadContents

        return ReadContents(self)

    def __getitem__(self, item):
        from janis_core.operators.operator import IndexOperator

        return IndexOperator(self, item)

    def basename(self):
        from .standard import BasenameOperator

        outtype = get_instantiated_type(self.returntype()).received_type()
        if not isinstance(outtype, (File, Directory)):
            raise Exception(
                "Basename operator can only be applied to steps of File / Directory type, received: "
                + str(outtype)
            )

        return BasenameOperator(self)

    def file_size(self):
        from .standard import FileSizeOperator

        return FileSizeOperator(self)


SelectorOrValue = Union[Selector, int, str, float]


class InputSelector(Selector):
    def __init__(self, input_to_select, remove_file_extension=None, **kwargs):
        # maybe worth validating the input_to_select identifier
        self.input_to_select = input_to_select

        if "use_basename" in kwargs:
            use_basename = kwargs["use_basename"]
            if remove_file_extension is None:
                remove_file_extension = use_basename
            Logger.warn(
                f"The 'use_basename' key is deprecated, please use 'remove_file_extension' instead: "
                f'InputSelector("{self.input_to_select}", remove_file_extension={str(use_basename)})'
            )

        self.remove_file_extension = remove_file_extension

    def returntype(self):
        # Todo: Work out how this can be achieved
        return File

    def to_string_formatter(self):
        kwarg = {self.input_to_select: self}
        from janis_core.operators.stringformatter import StringFormatter

        return StringFormatter(f"{{{self.input_to_select}}}", **kwarg)

    def __str__(self):
        return "inputs." + self.input_to_select

    def __repr__(self):
        return str(self)


class InputNodeSelector(Selector):
    def __init__(self, input_node):
        from janis_core.workflow.workflow import InputNode

        if input_node.node_type != NodeType.INPUT:  # input
            raise Exception(
                f"Error when creating InputOperator, '{input_node.id()}' was not an input node"
            )

        self.input_node: InputNode = input_node

    def id(self):
        return self.input_node.id()

    def returntype(self):
        out = first_value(self.input_node.outputs()).outtype

        if self.input_node is not None:
            import copy

            out = copy.copy(out)
            out.optional = False

        return out

    def __repr__(self):
        return "inputs." + self.input_node.id()

    def to_string_formatter(self):
        from janis_core.operators.stringformatter import StringFormatter

        key = self.input_node.id()
        kwarg = {key: self}
        return StringFormatter(f"{{{key}}}", **kwarg)


class StepOutputSelector(Selector):
    def __init__(self, node, tag):

        outputs = node.outputs()
        if tag not in outputs:
            raise TypeError(
                f"The step node {node.id()} did not have an output called '{tag}', "
                f"expected one of: {', '.join(outputs.keys())}"
            )

        self.node = node
        self.tag = tag

    def returntype(self):
        retval = self.node.outputs()[self.tag].outtype

        if hasattr(self.node, "scatter") and self.node.scatter:
            retval = Array(retval)
        return retval

    @staticmethod
    def from_tuple(step_tuple):
        return StepOutputSelector(step_tuple[0], step_tuple[1])

    def __repr__(self):
        parts = [p for p in (self.node.id(), self.tag) if p is not None]

        return ".".join(parts)

    def to_string_formatter(self):
        from janis_core.operators.stringformatter import StringFormatter

        key = self.node.id() + "_" + self.tag
        kwarg = {key: self}
        return StringFormatter(f"{{{key}}}", **kwarg)


class WildcardSelector(Selector):
    def __init__(self, wildcard, select_first=False):
        self.wildcard = wildcard
        self.select_first = select_first

    def returntype(self):
        return Array(Union[File, Directory])

    def to_string_formatter(self):
        raise Exception("A wildcard selector cannot be coerced into a StringFormatter")


class MemorySelector(InputSelector):
    def __init__(self):
        super().__init__("runtime_memory")

    def returntype(self):
        return Int(optional=True)

    def __str__(self):
        return "runtime_memory"

    def __repr__(self):
        return str(self)


class CpuSelector(InputSelector):
    def __init__(self, default=1):
        super().__init__("runtime_cpu")
        self.default = default

    def returntype(self):
        return Int(optional=bool(self.default is None))

    def __str__(self):
        return "runtime_cpu"

    def __repr__(self):
        return str(self)


class DiskSelector(InputSelector):
    def __init__(self):
        super().__init__("runtime_disks")

    def returntype(self):
        return Int(optional=True)

    def __str__(self):
        return "runtime_disks"

    def __repr__(self):
        return str(self)


class TimeSelector(InputSelector):
    def __init__(self, default=86400):
        """
        Specified in seconds
        :param default:
        """
        super().__init__("runtime_seconds")
        self.default = default

    def returntype(self):
        return Int(optional=bool(self.default is None))

    def __str__(self):
        return "runtime_seconds"

    def __repr__(self):
        return str(self)
