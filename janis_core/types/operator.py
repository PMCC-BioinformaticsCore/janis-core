from abc import ABC, abstractmethod
from typing import Union

from janis_core.types.selectors import Selector


class Operator(Selector, ABC):
    def __neg__(self):
        return NotOperator(self)

    def __and__(self, other):
        return AndOperator(self, other)

    def __rand__(self, other):
        return AndOperator(other, self)

    def __or__(self, other):
        return OrOperator(self, other)

    def __ror__(self, other):
        return OrOperator(other, self)

    def __add__(self, other):
        return AddOperator(self, other)

    def __radd__(self, other):
        return AddOperator(other, self)

    def __sub__(self, other):
        return SubtractOperator(self, other)

    def __rsub__(self, other):
        return SubtractOperator(other, self)

    def __mul__(self, other):
        return MultiplyOperator(self, other)

    def __rmul__(self, other):
        return MultiplyOperator(other, self)

    def __truediv__(self, other):
        return DivideOperator(self, other)

    def __rtruediv__(self, other):
        return DivideOperator(other, self)

    def __eq__(self, other):
        return EqualityOperator(self, other)

    def __ne__(self, other):
        return EqualityOperator(self, other)

    def __gt__(self, other):
        return GtOperator(self, other)

    def __ge__(self, other):
        return GteOperator(self, other)

    def __lt__(self, other):
        return LtOperator(self, other)

    def __le__(self, other):
        return LteOperator(self, other)

    def as_str(self):
        return AsStringOperator(self)

    def as_bool(self):
        return AsBoolOperator(self)

    def as_int(self):
        return AsIntOperator(self)

    def op_and(self, other):
        return AndOperator(self, other)

    def op_or(self, other):
        return OrOperator(self, other)


class InputOperator(Operator):
    def __init__(self, input_node):
        if input_node.node_type != 1:  # input
            raise Exception(
                f"Error when creating InputOperator, '{input_node.id()}' was not an input node"
            )

        self.input_node = input_node

    def __repr__(self):
        return "inputs." + self.input_node.id()


class StepOperator(Operator):
    def __init__(self, node, tag):
        self.node = node
        self.tag = tag

    @staticmethod
    def from_tuple(step_tuple):
        return StepOperator(step_tuple[0], step_tuple[1])

    def __repr__(self):
        return self.node.id() + "." + self.tag


OperatorOrValue = Union[Operator, Selector, int, str, float]


class SingleValueOperator(Operator, ABC):
    @staticmethod
    @abstractmethod
    def symbol():
        pass

    @staticmethod
    @abstractmethod
    def wdl_symbol():
        pass

    @staticmethod
    @abstractmethod
    def cwl_symbol():
        pass

    def __str__(self):
        return f"{self.symbol()}({self.internal})"

    def __init__(self, internal: OperatorOrValue):
        self.internal = internal


class TwoValueOperator(Operator, ABC):
    @staticmethod
    @abstractmethod
    def symbol():
        pass

    @staticmethod
    @abstractmethod
    def wdl_symbol():
        pass

    @staticmethod
    @abstractmethod
    def cwl_symbol():
        pass

    def __init__(self, lhs: OperatorOrValue, rhs: OperatorOrValue):
        self.lhs = lhs
        self.rhs = rhs

    def __str__(self):
        return f"({self.lhs} {self.symbol()} {self.rhs})"


# As operators


class AsStringOperator(SingleValueOperator):
    @staticmethod
    def symbol():
        return "str"

    @staticmethod
    def wdl_symbol():
        return "str"

    @staticmethod
    def cwl_symbol():
        return "str"


class AsBoolOperator(SingleValueOperator):
    @staticmethod
    def symbol():
        return "bool"

    @staticmethod
    def wdl_symbol():
        return "bool"

    @staticmethod
    def cwl_symbol():
        return "bool"


class AsIntOperator(SingleValueOperator):
    @staticmethod
    def symbol():
        return "int"

    @staticmethod
    def wdl_symbol():
        return "int"

    @staticmethod
    def cwl_symbol():
        return "int"


# Other single value operators


class NotOperator(SingleValueOperator):
    @staticmethod
    def symbol():
        return "!"

    @staticmethod
    def wdl_symbol():
        return "!"

    @staticmethod
    def cwl_symbol():
        return "!"


# Two value operators


class AndOperator(TwoValueOperator):
    @staticmethod
    def symbol():
        return "&&"

    @staticmethod
    def wdl_symbol():
        return "&&"

    @staticmethod
    def cwl_symbol():
        return "&&"


class OrOperator(TwoValueOperator):
    @staticmethod
    def symbol():
        return "||"

    @staticmethod
    def wdl_symbol():
        return "||"

    @staticmethod
    def cwl_symbol():
        return "||"


class EqualityOperator(TwoValueOperator):
    @staticmethod
    def symbol():
        return "=="

    @staticmethod
    def wdl_symbol():
        return "=="

    @staticmethod
    def cwl_symbol():
        return "=="


class InequalityOperator(TwoValueOperator):
    @staticmethod
    def symbol():
        return "!="

    @staticmethod
    def wdl_symbol():
        return "!="

    @staticmethod
    def cwl_symbol():
        return "!="


class GtOperator(TwoValueOperator):
    @staticmethod
    def symbol():
        return ">"

    @staticmethod
    def wdl_symbol():
        return ">"

    @staticmethod
    def cwl_symbol():
        return ">"


class GteOperator(TwoValueOperator):
    @staticmethod
    def symbol():
        return ">="

    @staticmethod
    def wdl_symbol():
        return ">="

    @staticmethod
    def cwl_symbol():
        return ">="


class LtOperator(TwoValueOperator):
    @staticmethod
    def symbol():
        return "<"

    @staticmethod
    def wdl_symbol():
        return "<"

    @staticmethod
    def cwl_symbol():
        return "<"


class LteOperator(TwoValueOperator):
    @staticmethod
    def symbol():
        return "<="

    @staticmethod
    def wdl_symbol():
        return "<="

    @staticmethod
    def cwl_symbol():
        return "<="


class AddOperator(TwoValueOperator):
    @staticmethod
    def symbol():
        return "+"

    @staticmethod
    def wdl_symbol():
        return "+"

    @staticmethod
    def cwl_symbol():
        return "+"


class SubtractOperator(TwoValueOperator):
    @staticmethod
    def symbol():
        return "-"

    @staticmethod
    def wdl_symbol():
        return "-"

    @staticmethod
    def cwl_symbol():
        return "-"


class MultiplyOperator(TwoValueOperator):
    @staticmethod
    def symbol():
        return "*"

    @staticmethod
    def wdl_symbol():
        return "*"

    @staticmethod
    def cwl_symbol():
        return "*"


class DivideOperator(TwoValueOperator):
    @staticmethod
    def symbol():
        return "/"

    @staticmethod
    def wdl_symbol():
        return "/"

    @staticmethod
    def cwl_symbol():
        return "/"
