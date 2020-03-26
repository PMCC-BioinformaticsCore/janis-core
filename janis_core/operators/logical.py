from abc import ABC, abstractmethod
from typing import List
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
