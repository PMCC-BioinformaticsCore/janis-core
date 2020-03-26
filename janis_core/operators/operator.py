from abc import ABC, abstractmethod
from typing import Union, List

from janis_core.types import DataType, get_instantiated_type

from janis_core.utils import first_value

from janis_core.operators.selectors import Selector


class Operator(Selector, ABC):
    def __neg__(self):
        from janis_core.operators.logical import NotOperator

        return NotOperator(self)

    def __and__(self, other):
        from janis_core.operators.logical import AndOperator

        return AndOperator(self, other)

    def __rand__(self, other):
        from janis_core.operators.logical import AndOperator

        return AndOperator(other, self)

    def __or__(self, other):
        from janis_core.operators.logical import OrOperator

        return OrOperator(self, other)

    def __ror__(self, other):
        from janis_core.operators.logical import OrOperator

        return OrOperator(other, self)

    def __add__(self, other):
        from janis_core.operators.logical import AddOperator

        return AddOperator(self, other)

    def __radd__(self, other):
        from janis_core.operators.logical import AddOperator

        return AddOperator(other, self)

    def __sub__(self, other):
        from janis_core.operators.logical import SubtractOperator

        return SubtractOperator(self, other)

    def __rsub__(self, other):
        from janis_core.operators.logical import SubtractOperator

        return SubtractOperator(other, self)

    def __mul__(self, other):
        from janis_core.operators.logical import MultiplyOperator

        return MultiplyOperator(self, other)

    def __rmul__(self, other):
        from janis_core.operators.logical import MultiplyOperator

        return MultiplyOperator(other, self)

    def __truediv__(self, other):
        from janis_core.operators.logical import DivideOperator

        return DivideOperator(self, other)

    def __rtruediv__(self, other):
        from janis_core.operators.logical import DivideOperator

        return DivideOperator(other, self)

    def __eq__(self, other):
        from janis_core.operators.logical import EqualityOperator

        return EqualityOperator(self, other)

    def __ne__(self, other):
        from janis_core.operators.logical import EqualityOperator

        return EqualityOperator(self, other)

    def __gt__(self, other):
        from janis_core.operators.logical import GtOperator

        return GtOperator(self, other)

    def __ge__(self, other):
        from janis_core.operators.logical import GteOperator

        return GteOperator(self, other)

    def __lt__(self, other):
        from janis_core.operators.logical import LtOperator

        return LtOperator(self, other)

    def __le__(self, other):
        from janis_core.operators.logical import LteOperator

        return LteOperator(self, other)

    def as_str(self):
        return AsStringOperator(self)

    def as_bool(self):
        return AsBoolOperator(self)

    def as_int(self):
        return AsIntOperator(self)

    def op_and(self, other):
        from janis_core.operators.logical import AndOperator

        return AndOperator(self, other)

    def op_or(self, other):
        from janis_core.operators.logical import OrOperator

        return OrOperator(self, other)


class InputOperator(Operator):
    def __init__(self, input_node):
        if input_node.node_type != 1:  # input
            raise Exception(
                f"Error when creating InputOperator, '{input_node.id()}' was not an input node"
            )

        self.input_node = input_node

    def returntype(self):
        return first_value(self.input_node.inputs()).intype

    def __repr__(self):
        return "inputs." + self.input_node.id()


class StepOperator(Operator):
    def __init__(self, node, tag):
        self.node = node
        self.tag = tag

    def returntype(self):
        return self.node.inputs()[self.tag].intype

    @staticmethod
    def from_tuple(step_tuple):
        return StepOperator(step_tuple[0], step_tuple[1])

    def __repr__(self):
        return self.node.id() + "." + self.tag

    def as_operator(self):
        return self


OperatorOrValue = Union[Operator, Selector, int, str, float]


class FunctionOperator(Operator, ABC):
    def __init__(self, *args):
        self.args: List[Selector] = args

    @abstractmethod
    def argtypes(self) -> List[DataType]:
        pass

    def validate(self):
        args = self.args
        argtypes = self.argtypes()

        if len(args) < len(argtypes):
            # missing params
            nmissing = len(argtypes) - len(args)
            raise TypeError(
                f"{self.__class__.__name__} missing {nmissing} required positional argument"
            )
        elif len(args) > len(argtypes):
            raise TypeError(
                f"{self.__class__.__name__} expects {len(argtypes)} arguments, but {len(args)} were given"
            )

        errors = []
        for i in range(len(args)):
            expected = get_instantiated_type(argtypes[i])
            received = get_instantiated_type(args[i])
            if not expected.can_receive_from(received):
                errors.append(
                    f"argument {i+1} '{received.id()}' was incompatible with the expected {expected.id()}"
                )

        if errors:
            singularprefix = "were errors" if len(errors) != 1 else "was an error"
            raise TypeError(
                f"There {singularprefix} in the arguments when validating '{self.__class__.__name__}', "
                + ", ".join(errors)
            )

        return True

    @abstractmethod
    def to_wdl(self, unwrap_operator, *args):
        pass

    @abstractmethod
    def to_cwl(self, unwrap_operator, *args):
        pass


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
