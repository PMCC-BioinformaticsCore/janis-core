from abc import ABC, abstractmethod
from typing import List, Union
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
from .operator import Operator, OperatorOrValue, SingleValueOperator, TwoValueOperator


# As operators


def or_prev_conds(prevconditions: List[Operator]):
    if len(prevconditions) == 0:
        return None
    elif len(prevconditions) == 1:
        return prevconditions[0]
    elif len(prevconditions) == 2:
        return OrOperator(prevconditions[0], prevconditions[1])
    return OrOperator(prevconditions[0], or_prev_conds(prevconditions[1:]))


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
        lhs: DataType = self.args[0].returntype()
        rhs: DataType = self.args[0].returntype()

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
