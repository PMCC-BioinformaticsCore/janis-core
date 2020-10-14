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
    @staticmethod
    def friendly_signature():
        return "Any? -> Boolean"

    def argtypes(self) -> List[DataType]:
        return [AnyType]

    def returntype(self):
        return Boolean

    def __str__(self):
        arg = self.args[0]
        return f"isdefined({arg})"

    def __repr__(self):
        return str(self)

    def evaluate(self, inputs):
        return self.evaluate_arg(self.args[0], inputs) is not None

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

    @staticmethod
    def friendly_signature():
        return "(condition: Boolean, value_if_true: X, value_if_false: Y) -> Union[X,Y]"

    def argtypes(self) -> List[DataType]:
        return [Boolean, AnyType, AnyType]

    def returntype(self):
        args = []
        for a in self.args[1:]:
            if isinstance(a, Selector):
                args.append(get_instantiated_type(a.returntype()))
            else:
                args.append(get_instantiated_type(a))
        return UnionType(*args)

    def __str__(self):
        cond, v1, v2 = self.args
        return f"{cond} ? {v1} : {v2}"

    def __repr__(self):
        return str(self)

    def evaluate(self, inputs):
        cond, iftrue, iffalse = self.args
        result = iftrue if self.evaluate_arg(cond, inputs) else iffalse
        return self.evaluate_arg(result, inputs)

    def to_wdl(self, unwrap_operator, *args):
        cond, v1, v2 = [unwrap_operator(a) for a in self.args]
        return f"if ({cond}) then {v1} else {v2}"

    def to_cwl(self, unwrap_operator, *args):
        cond, v1, v2 = [unwrap_operator(a) for a in self.args]
        return f"{cond} ? {v1} : {v2}"


class AssertNotNull(Operator):
    @staticmethod
    def friendly_signature():
        return "X? -> X"

    def argtypes(self) -> List[DataType]:
        return [AnyType]

    def evaluate(self, inputs):
        result = self.evaluate_arg(self.args[0], inputs)
        assert result is not None
        return result

    def to_wdl(self, unwrap_operator, *args):
        arg = unwrap_operator(self.args[0])
        return f"select_first([{arg}])"

    def to_cwl(self, unwrap_operator, *args):
        return unwrap_operator(self.args[0])

    def returntype(self):
        from copy import copy

        a = self.args[0]
        if isinstance(a, Selector):
            ret = get_instantiated_type(a.returntype())
        else:
            ret = get_instantiated_type(a)

        ret = copy(ret)
        ret.optional = False
        return ret


# Other single value operators


class NotOperator(SingleValueOperator):
    @staticmethod
    def friendly_signature():
        return "Boolean -> Boolean"

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

    @staticmethod
    def apply_to(value):
        return not value

    # Two value operators


class AndOperator(TwoValueOperator):
    @staticmethod
    def friendly_signature():
        return "Boolean, Boolean -> Boolean"

    @staticmethod
    def symbol():
        return "&&"

    @staticmethod
    def wdl_symbol():
        return "&&"

    @staticmethod
    def cwl_symbol():
        return "&&"

    @staticmethod
    def apply_to(arg1, arg2):
        return arg1 and arg2

    def returntype(self):
        return Boolean

    def argtypes(self):
        return [Boolean, Boolean]


class OrOperator(TwoValueOperator):
    @staticmethod
    def friendly_signature():
        return "Boolean, Boolean -> Boolean"

    @staticmethod
    def symbol():
        return "||"

    @staticmethod
    def wdl_symbol():
        return "||"

    @staticmethod
    def cwl_symbol():
        return "||"

    @staticmethod
    def apply_to(arg1, arg2):
        return arg1 or arg2

    def returntype(self):
        return Boolean

    def argtypes(self):
        return [Boolean, Boolean]


class EqualityOperator(TwoValueOperator):
    @staticmethod
    def friendly_signature():
        return "Any, Any -> Boolean"

    @staticmethod
    def symbol():
        return "=="

    @staticmethod
    def wdl_symbol():
        return "=="

    @staticmethod
    def cwl_symbol():
        return "=="

    @staticmethod
    def apply_to(arg1, arg2):
        return arg1 == arg2

    def returntype(self):
        return Boolean

    def argtypes(self):
        return [AnyType, AnyType]


class InequalityOperator(TwoValueOperator):
    @staticmethod
    def friendly_signature():
        return "Any, Any -> Boolean"

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

    @staticmethod
    def apply_to(arg1, arg2):
        return arg1 != arg2


class GtOperator(TwoValueOperator):
    @staticmethod
    def friendly_signature():
        return "Comparable, Comparable -> Boolean"

    @staticmethod
    def symbol():
        return ">"

    @staticmethod
    def wdl_symbol():
        return ">"

    @staticmethod
    def cwl_symbol():
        return ">"

    @staticmethod
    def apply_to(arg1, arg2):
        return arg1 > arg2

    def returntype(self):
        return Boolean

    def argtypes(self):
        return [Union[String, Int, Float, Double], Union[String, Int, Float, Double]]


class GteOperator(TwoValueOperator):
    @staticmethod
    def friendly_signature():
        return "Comparable, Comparable -> Boolean"

    @staticmethod
    def symbol():
        return ">="

    @staticmethod
    def wdl_symbol():
        return ">="

    @staticmethod
    def cwl_symbol():
        return ">="

    @staticmethod
    def apply_to(arg1, arg2):
        return arg1 >= arg2

    def returntype(self):
        return Boolean

    def argtypes(self):
        return [Union[String, Int, Float, Double], Union[String, Int, Float, Double]]


class LtOperator(TwoValueOperator):
    @staticmethod
    def friendly_signature():
        return "Comparable, Comparable -> Boolean"

    @staticmethod
    def symbol():
        return "<"

    @staticmethod
    def wdl_symbol():
        return "<"

    @staticmethod
    def cwl_symbol():
        return "<"

    @staticmethod
    def apply_to(arg1, arg2):
        return arg1 < arg2

    def returntype(self):
        return Boolean

    def argtypes(self):
        return [Union[String, Int, Float, Double], Union[String, Int, Float, Double]]


class LteOperator(TwoValueOperator):
    @staticmethod
    def friendly_signature():
        return "Comparable, Comparable -> Boolean"

    @staticmethod
    def symbol():
        return "<="

    @staticmethod
    def wdl_symbol():
        return "<="

    @staticmethod
    def cwl_symbol():
        return "<="

    @staticmethod
    def apply_to(arg1, arg2):
        return arg1 <= arg2

    def returntype(self):
        return Boolean

    def argtypes(self):
        return [Union[String, Int, Float, Double], Union[String, Int, Float, Double]]


class AddOperator(TwoValueOperator):
    @staticmethod
    def friendly_signature():
        return "Any, Any -> Any"

    @staticmethod
    def symbol():
        return "+"

    @staticmethod
    def wdl_symbol():
        return "+"

    @staticmethod
    def cwl_symbol():
        return "+"

    @staticmethod
    def apply_to(arg1, arg2):
        return arg1 + arg2

    def argtypes(self):
        return [AnyType, AnyType]

    def returntype(self):
        lhs_val: DataType = self.args[0]
        rhs_val: DataType = self.args[1]

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
    def friendly_signature():
        return "NumericType, NumericType -> NumericType"

    @staticmethod
    def symbol():
        return "-"

    @staticmethod
    def wdl_symbol():
        return "-"

    @staticmethod
    def cwl_symbol():
        return "-"

    @staticmethod
    def apply_to(arg1, arg2):
        return arg1 - arg2

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
    def friendly_signature():
        return "Numeric, NumericType -> NumericType"

    @staticmethod
    def symbol():
        return "*"

    @staticmethod
    def wdl_symbol():
        return "*"

    @staticmethod
    def cwl_symbol():
        return "*"

    @staticmethod
    def apply_to(arg1, arg2):
        return arg1 * arg2

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
    def friendly_signature():
        return "Numeric, NumericType -> NumericType"

    @staticmethod
    def symbol():
        return "/"

    @staticmethod
    def wdl_symbol():
        return "/"

    @staticmethod
    def cwl_symbol():
        return "/"

    @staticmethod
    def apply_to(arg1, arg2):
        return arg1 / arg2

    def argtypes(self):
        return [NumericType, NumericType]

    def returntype(self):
        if any(isinstance(t, Double) for t in self.args):
            return Double
        return Float


class FloorOperator(Operator):
    @staticmethod
    def friendly_signature():
        return "Numeric, NumericType -> Int"

    def argtypes(self) -> List[DataType]:
        return [NumericType]

    def returntype(self):
        return Int

    def __str__(self):
        arg = self.args[0]
        return f"floor({arg})"

    def __repr__(self):
        return str(self)

    def to_wdl(self, unwrap_operator, *args):
        arg = unwrap_operator(self.args[0])
        return f"floor({arg})"

    def to_cwl(self, unwrap_operator, *args):
        arg = unwrap_operator(self.args[0])
        return f"Math.floor({arg})"

    def evaluate(self, inputs):
        from math import floor

        result = self.evaluate_arg(self.args[0], inputs)
        return floor(result)


class CeilOperator(Operator):
    @staticmethod
    def friendly_signature():
        return "Numeric, NumericType -> Int"

    def argtypes(self) -> List[DataType]:
        return [NumericType]

    def returntype(self):
        return Int

    def __str__(self):
        arg = self.args[0]
        return f"ceil({arg})"

    def __repr__(self):
        return str(self)

    def to_wdl(self, unwrap_operator, *args):
        arg = unwrap_operator(self.args[0])
        return f"ceil({arg})"

    def to_cwl(self, unwrap_operator, *args):
        arg = unwrap_operator(self.args[0])
        return f"Math.ceil({arg})"

    def evaluate(self, inputs):
        from math import ceil

        result = self.evaluate_arg(self.args[0], inputs)
        return ceil(result)


class RoundOperator(Operator):
    @staticmethod
    def friendly_signature():
        return "Numeric, NumericType -> Int"

    def argtypes(self) -> List[DataType]:
        return [NumericType]

    def returntype(self):
        return Int

    def __str__(self):
        arg = self.args[0]
        return f"round({arg})"

    def __repr__(self):
        return str(self)

    def to_wdl(self, unwrap_operator, *args):
        arg = unwrap_operator(self.args[0])
        return f"round({arg})"

    def to_cwl(self, unwrap_operator, *args):
        arg = unwrap_operator(self.args[0])
        return f"Math.round({arg})"

    def evaluate(self, inputs):
        result = self.evaluate_arg(self.args[0], inputs)
        return round(result)
