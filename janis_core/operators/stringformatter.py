from typing import Optional, List, Dict, Tuple

from janis_core.utils import first_value

from janis_core.types import String, AnyType
from janis_core.operators.logical import Operator, AddOperator
from janis_core.utils.bracketmatching import get_keywords_between_braces
from janis_core.utils.errors import (
    TooManyArgsException,
    IncorrectArgsException,
    InvalidByProductException,
    ConflictingArgumentsException,
)
from janis_core.utils.logger import Logger


class StringFormatter(Operator):
    def returntype(self):
        return String()

    def argtypes(self):
        return [String, Optional[AnyType]]

    @staticmethod
    def friendly_signature():
        return "String, **kwargs -> String"

    def validate(self, perform_typecheck=False):
        return True

    def __init__(self, format: str, **kwargs):
        super().__init__([])
        # ignore super().__init__ call
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

    def to_cwl(self, unwrap_operator, *args):
        raise Exception("Don't use this method")

    def to_wdl(self, unwrap_operator, *args):
        raise Exception("Don't use this method")

    def evaluate(self, inputs):
        resolvedvalues = {
            k: self.evaluate_arg(v, inputs) for k, v in self.kwargs.items()
        }

        values_that_are_lists = {
            k: v for k, v in resolvedvalues.items() if isinstance(v, list)
        }

        inp_combinations: List[dict] = [{}]

        if len(values_that_are_lists) > 0:
            l = len(first_value(values_that_are_lists))
            list_values_that_are_different = sum(
                0 if len(v) == l else 1 for v in values_that_are_lists.values()
            )

            if list_values_that_are_different == 0:
                # dot product
                inp_combinations = [
                    {k: v[i] for k, v in values_that_are_lists.items()}
                    for i in range(l)
                ]
            elif list_values_that_are_different == 1:
                # cross product
                inp_combinations = self.generate_combinations_of_input_dicts(
                    values_that_are_lists=list(values_that_are_lists.items())
                )
            else:
                l_lengths = ", ".join(
                    f"{k}={len(v)}" for k, v in values_that_are_lists.items()
                )
                raise Exception(
                    "String Formatter evaluation doesn't support scattering for list of "
                )

        evaluated_combinations = [
            self.resolve_with_resolved_values(**{**resolvedvalues, **c})
            for c in inp_combinations
        ]
        if len(evaluated_combinations) == 0:
            raise Exception(
                "Something happened when resolving inputs with input values "
                + str(inputs)
            )
        elif len(evaluated_combinations) == 1:
            return evaluated_combinations[0]
        else:
            return evaluated_combinations

    def rewrite_operator(self, args_to_rewrite: dict):
        return self.__class__(
            self._format, **self.substitute_arg(args_to_rewrite, self.kwargs)
        )

    @staticmethod
    def generate_combinations_of_input_dicts(
        values_that_are_lists: List[Tuple[str, List[any]]]
    ) -> List[Dict]:

        if len(values_that_are_lists) == 0:
            return []
        key = values_that_are_lists[0][0]
        values = values_that_are_lists[0][1]

        if len(values_that_are_lists) == 1:
            return [{key: v} for v in values]

        combinations = []
        for v in values:
            for c in StringFormatter.generate_combinations_of_input_dicts(
                values_that_are_lists[1:]
            ):
                combinations.append({**c, key: v})

        return combinations

    def __repr__(self):
        val = self._format
        for k, v in self.kwargs.items():
            val = val.replace(f"{{{k}}}", f"{{{str(v)}}}")
        return val

    def get_leaves(self):
        leaves = []
        for a in self.kwargs.values():
            if isinstance(a, Operator):
                leaves.extend(a.get_leaves())
            else:
                leaves.append(a)
        return leaves

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
        from janis_core.operators.selectors import InputSelector

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

    def to_string_formatter(self):
        return self
