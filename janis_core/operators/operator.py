from abc import ABC, abstractmethod
from typing import List, Union

from janis_core.operators.selectors import Selector, StringFormatter
from janis_core.types import DataType, get_instantiated_type, UnionType
from janis_core.types.common_data_types import String, Boolean, Int, AnyType, Array


class Operator(Selector, ABC):
    def __init__(self, *args):
        self.args: List[Union[Selector, any]] = list(args)

    def get_leaves(self):
        leaves = []
        for a in self.args:
            if isinstance(a, Operator):
                leaves.extend(a.get_leaves())
            else:
                leaves.append(a)
        return leaves

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
            received_arg = args[i]
            if isinstance(received_arg, Selector):
                received = get_instantiated_type(args[i].returntype())
            else:
                received = get_instantiated_type(received_arg)
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

    def to_string_formatter(self):
        import re

        key = re.sub(r"\W+", "", str(self))
        kwarg = {key: self}
        return StringFormatter(f"{{{key}}}", **kwarg)


class IndexOperator(Operator, ABC):
    def __init__(self, base, index):
        super().__init__(base, index)

    def argtypes(self):
        return [Array(AnyType), Int]

    def returntype(self):
        return self.args[0].returntype().subtype()

    def __str__(self):
        base, index = self.args
        return f"{base}[{index}]"

    def __repr__(self):
        return str(self)

    def to_wdl(self, unwrap_operator, *args):
        base, index = self.args
        return f"{base}[{index}]"

    def to_cwl(self, unwrap_operator, *args):
        base, index = self.args
        return f"{base}[{index}]"


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
        internal = self.args[0]
        return f"{self.symbol()}({internal})"

    def to_wdl(self, unwrap_operator, *args):
        return f"{self.symbol()}({unwrap_operator(*args)})"

    def to_cwl(self, unwrap_operator, *args):
        return f"{self.symbol()}({unwrap_operator(*args)})"


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

    def to_wdl(self, unwrap_operator, *args):
        arg1, arg2 = [unwrap_operator(a) for a in self.args]
        return f"({arg1} {self.wdl_symbol()} {arg2})"

    def to_cwl(self, unwrap_operator, *args):
        arg1, arg2 = [unwrap_operator(a) for a in self.args]
        return f"({arg1} {self.cwl_symbol()} {arg2})"

    def __str__(self):
        args = self.args
        return f"({args[0]} {self.symbol()} {args[1]})"


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

    def argtypes(self):
        return [AnyType]

    def returntype(self):
        return String


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

    def argtypes(self):
        return [AnyType]

    def returntype(self):
        return Boolean


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

    def argtypes(self):
        return [AnyType]

    def returntype(self):
        return Int
