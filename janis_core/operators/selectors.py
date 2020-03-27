from typing import Union

from janis_core.utils import first_value

from janis_core.utils.errors import (
    TooManyArgsException,
    IncorrectArgsException,
    InvalidByProductException,
    ConflictingArgumentsException,
)
from janis_core.utils.logger import Logger
from janis_core.types.common_data_types import Array, String, File, Directory, Int

from janis_core.utils.bracketmatching import get_keywords_between_braces
from abc import ABC, abstractmethod


class Selector(ABC):
    @abstractmethod
    def returntype(self):
        pass

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
        from janis_core.operators.operator import AsStringOperator

        return AsStringOperator(self)

    def as_bool(self):
        from janis_core.operators.operator import AsBoolOperator

        return AsBoolOperator(self)

    def as_int(self):
        from janis_core.operators.operator import AsIntOperator

        return AsIntOperator(self)

    def op_and(self, other):
        from janis_core.operators.logical import AndOperator

        return AndOperator(self, other)

    def op_or(self, other):
        from janis_core.operators.logical import OrOperator

        return OrOperator(self, other)

    def __getitem__(self, item):
        from janis_core.operators.operator import IndexOperator

        return IndexOperator(self, item)


class InputSelector(Selector):
    def __init__(self, input_to_select, use_basename=None):
        # maybe worth validating the input_to_select identifier
        self.input_to_select = input_to_select
        self.use_basename = use_basename

    def returntype(self):
        # Todo: Work out how this can be achieved
        return File

    def to_string_formatter(self):
        kwarg = {self.input_to_select: self}
        return StringFormatter(f"{{{self.input_to_select}}}", **kwarg)

    def __radd__(self, other):
        return StringFormatter(other) + self.to_string_formatter()

    def __add__(self, other):
        return self.to_string_formatter() + other


class InputNodeSelector(Selector):
    def __init__(self, input_node):
        if input_node.node_type != 1:  # input
            raise Exception(
                f"Error when creating InputOperator, '{input_node.id()}' was not an input node"
            )

        self.input_node = input_node

    def returntype(self):
        out = first_value(self.input_node.outputs()).outtype

        if self.input_node is not None:
            import copy

            out = copy.copy(out)
            out.optional = False

        return out

    def __repr__(self):
        return "inputs." + self.input_node.id()


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
        if self.node.scatter:
            retval = Array(retval)
        return retval

    @staticmethod
    def from_tuple(step_tuple):
        return StepOutputSelector(step_tuple[0], step_tuple[1])

    def __repr__(self):
        return self.node.id() + "." + self.tag

    def as_operator(self):
        return self


class WildcardSelector(Selector):
    def __init__(self, wildcard):
        self.wildcard = wildcard

    def returntype(self):
        return Array(Union[File, Directory])


class MemorySelector(InputSelector):
    def __init__(self):
        super().__init__("runtime_memory")

    def returntype(self):
        return Int(optional=True)


class CpuSelector(InputSelector):
    def __init__(self, default=1):
        super().__init__("runtime_cpu")
        self.default = default

    def returntype(self):
        return Int(optional=bool(self.default is None))


class StringFormatter(Selector):
    def returntype(self):
        return String()

    def __init__(self, format: str, **kwargs):
        self._format: str = format

        keywords, balance = get_keywords_between_braces(self._format)

        if balance > 0:
            Logger.warn(
                "There was an imbalance of braces in the string _format, this might cause issues with concatenation"
            )

        skwargs = set(kwargs.keys())

        if not keywords == skwargs:
            # what's the differences
            if not keywords.issubset(skwargs):
                raise IncorrectArgsException(
                    "The _format required additional arguments to be provided by "
                    "**kwargs, requires the keys:" + ", ".join(keywords - skwargs)
                )
            else:
                raise TooManyArgsException(
                    "The **kwargs contained unrecognised keys: "
                    + ", ".join(skwargs - keywords)
                )

        self.kwargs = kwargs

    resolved_types = [str, int, float]

    def resolve_with_resolved_values(self, **resolved_values):

        s1 = set(self.kwargs.keys())
        actual_keys, _ = get_keywords_between_braces(self._format)
        if s1 != actual_keys:
            diff = (actual_keys - s1).union(s1 - actual_keys)

            raise Exception(
                "The format for the string builder has changed since runtime, or an internal error has"
                " occurred. The following keys did not appear in both sets: "
                + ", ".join(diff)
            )

        s2 = set(resolved_values.keys())

        missing_keys = s1 - s2
        if len(missing_keys) > 0:
            raise IncorrectArgsException(
                "There were missing parameters when formatting string: "
                + ", ".join(missing_keys)
            )

        unresolved_values = [
            f"{r} ({type(resolved_values[r]).__name__})"
            for r in resolved_values
            if not any(
                isinstance(resolved_values[r], t)
                for t in StringFormatter.resolved_types
            )
        ]
        if len(unresolved_values) > 0:
            raise ValueError(
                "There were unresolved parameters when formatting string: "
                + ", ".join(unresolved_values)
            )

        retval = self._format
        for k in resolved_values:
            retval = retval.replace(f"{{{k}}}", str(resolved_values[k]))
        return retval

    def __radd__(self, other):
        return StringFormatter(other) + self

    def __add__(self, other):
        if isinstance(other, str):
            # check if it has args in it
            keywords = get_keywords_between_braces(other)
            if len(keywords) > 0:
                invalidkwargs = [k for k in self.kwargs if k not in self.kwargs]
                if len(invalidkwargs) > 0:
                    raise InvalidByProductException(
                        f"The string to be concatenated contained placeholder(s) ({', '.join(invalidkwargs)})"
                        f"that were not in the original StringFormatter"
                    )
            return self._create_new_formatter_from_strings_and_args(
                [self._format, other], **self.kwargs
            )

        elif isinstance(other, InputSelector):
            return self + other.to_string_formatter()

        elif isinstance(other, StringFormatter):
            # check if args overlap and they're different
            s1 = set(self.kwargs.keys())
            s2 = set(other.kwargs.keys())
            intersection = s1.intersection(s2)

            if len(intersection) > 0:
                not_same_args = [
                    k for k in intersection if self.kwargs[k] != other.kwargs[k]
                ]
                if len(not_same_args) > 0:
                    raise ConflictingArgumentsException(
                        f"Couldn't concatenate formats as there keys ({', '.join(not_same_args)}) "
                        f"that were not equal between formatters "
                    )

            # yeah we sweet
            new_args = {**self.kwargs, **other.kwargs}
            return StringFormatter._create_new_formatter_from_strings_and_args(
                [self._format, other._format], **new_args
            )

    @staticmethod
    def _create_new_formatter_from_strings_and_args(strings: [str], **kwargs):
        new_format = "".join(strings)
        try:
            return StringFormatter(new_format, **kwargs)
        except IncorrectArgsException as e:
            new_params = set(
                get_keywords_between_braces(new_format)[0] - set(kwargs.keys())
            )
            raise InvalidByProductException(
                "Joining the input files (to '{new_format}') created the new params: "
                + ", ".join(new_params)
            )
