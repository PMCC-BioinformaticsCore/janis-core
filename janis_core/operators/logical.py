from abc import ABC, abstractmethod
from typing import List, Union

from janis_core.types import get_instantiated_type

from ..types import UnionType
from ..types.common_data_types import (
    Boolean,
    String,
    Float,
    DataType,
    Int,
    File,
    Directory,
    Double,
    AnyType,
    NumericType,
)
from .operator import Selector, Operator, SingleValueOperator, TwoValueOperator


# As operators


def or_prev_conds(prevconditions: List[Operator]):
    if len(prevconditions) == 0:
        return None
    elif len(prevconditions) == 1:
        return prevconditions[0]
    elif len(prevconditions) == 2:
        return OrOperator(prevconditions[0], prevconditions[1])
    return OrOperator(prevconditions[0], or_prev_conds(prevconditions[1:]))


class IsDefined(Operator, ABC):
    def argtypes(self) -> List[DataType]:
        return [AnyType]

    def returntype(self):
        return Boolean

    def __str__(self):
        arg = self.args[0]
        return f"isdefined({arg})"

    def __repr__(self):
        return str(self)

    def to_cwl(self, unwrap_operator, *args):
        arg = unwrap_operator(self.args[0])
        # 2 equals (!=) in javascript will coerce undefined to equal null
        return f"({arg} != null)"

    def to_wdl(self, unwrap_operator, *args):
        arg = unwrap_operator(self.args[0])
        return f"defined({arg})"


class If(Operator, ABC):
    def __init__(self, condition, value_if_true, value_if_false):
        super().__init__(condition, value_if_true, value_if_false)

    def argtypes(self) -> List[DataType]:
        return [Boolean, AnyType, AnyType]

    def returntype(self):
        return UnionType(*self.args[1:])

    def __str__(self):
        cond, v1, v2 = self.args
        return f"{cond} ? {v1} : {v2}"

    def __repr__(self):
        return str(self)

    def to_wdl(self, unwrap_operator, *args):
        cond, v1, v2 = [unwrap_operator(a) for a in self.args]
        return f"if ({cond}) then {v1} else {v2}"

    def to_cwl(self, unwrap_operator, *args):
        cond, v1, v2 = [unwrap_operator(a) for a in self.args]
        return f"{cond} ? {v1} : {v2}"


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

    def argtypes(self):
        return [Boolean]

    def returntype(self):
        return Boolean


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

    def returntype(self):
        return Boolean

    def argtypes(self):
        return [Boolean, Boolean]


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

    def returntype(self):
        return Boolean

    def argtypes(self):
        return [Boolean, Boolean]


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

    def returntype(self):
        return Boolean

    def argtypes(self):
        return [AnyType, AnyType]


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

    def returntype(self):
        return Boolean

    def argtypes(self):
        return [AnyType, AnyType]


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

    def returntype(self):
        return Boolean

    def argtypes(self):
        return [Union[String, Int, Float, Double], Union[String, Int, Float, Double]]


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

    def returntype(self):
        return Boolean

    def argtypes(self):
        return [Union[String, Int, Float, Double], Union[String, Int, Float, Double]]


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

    def returntype(self):
        return Boolean

    def argtypes(self):
        return [Union[String, Int, Float, Double], Union[String, Int, Float, Double]]


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

    def returntype(self):
        return Boolean

    def argtypes(self):
        return [Union[String, Int, Float, Double], Union[String, Int, Float, Double]]


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

    def argtypes(self):
        return [AnyType, AnyType]

    def returntype(self):
        lhs_val: DataType = self.args[0]
        rhs_val: DataType = self.args[0]

        if isinstance(lhs_val, Selector):
            lhs = get_instantiated_type(lhs_val.returntype())
        else:
            lhs = get_instantiated_type(lhs_val)

        if isinstance(rhs_val, Selector):
            rhs = get_instantiated_type(rhs_val.returntype())
        else:
            rhs = get_instantiated_type(rhs_val)

        if isinstance(lhs, (String, File, Directory)) or isinstance(
            rhs, (String, File, Directory)
        ):
            return String
        if isinstance(lhs, Float) or isinstance(rhs, Float):
            return Double
        if isinstance(lhs, Float) or isinstance(rhs, Float):
            return Float
        if isinstance(lhs, Int) and isinstance(rhs, Int):
            return Int

        raise TypeError(f"Unsure how to derive returntype from {lhs.id()} + {rhs.id()}")


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

    def argtypes(self):
        return [Union[Int, Double, Float], Union[Int, Double, Float]]

    def returntype(self):
        if any(isinstance(t, Double) for t in self.args):
            return Double
        if any(isinstance(t, Float) for t in self.args):
            return Float
        return Int


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

    def argtypes(self):
        return [NumericType, NumericType]

    def returntype(self):
        if any(isinstance(t, Double) for t in self.args):
            return Double
        if any(isinstance(t, Float) for t in self.args):
            return Float
        return Int


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

    def argtypes(self):
        return [NumericType, NumericType]

    def returntype(self):
        if any(isinstance(t, Double) for t in self.args):
            return Double
        return Float
