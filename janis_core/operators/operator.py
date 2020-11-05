from abc import ABC, abstractmethod
from typing import List, Union

from janis_core.operators.selectors import Selector, InputSelector, InputNodeSelector
from janis_core.types import DataType, get_instantiated_type, Float
from janis_core.types.common_data_types import String, Boolean, Int, AnyType, Array


class Operator(Selector, ABC):
    def __init__(self, *args):
        self.args: List[Union[Selector, any]] = list(args)

        self.validate()

    def requires_contents(self):
        """
        A subclass should set this to TRUE
        :return:
        """
        return self._requires_contents(self.args)

    @staticmethod
    def _requires_contents(obj):

        if isinstance(obj, list):
            return any(Operator._requires_contents(o) for o in obj)
        elif isinstance(obj, Selector):
            return obj.requires_contents()
        return False

    @staticmethod
    @abstractmethod
    def friendly_signature():
        pass

    def get_leaves(self):
        leaves = []
        for a in self.args:
            if isinstance(a, Operator):
                leaves.extend(a.get_leaves())
            else:
                (leaves.extend if isinstance(a, list) else leaves.append)(a)
        return leaves

    @abstractmethod
    def argtypes(self) -> List[DataType]:
        pass

    def validate(self, perform_typecheck=False):
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

        if perform_typecheck:
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

    @staticmethod
    def evaluate_arg(arg, inputs):
        if arg is None:
            return None
        elif isinstance(arg, list):
            return [Operator.evaluate_arg(a, inputs) for a in arg]
        elif isinstance(arg, (str, int, float, bool)):
            return arg

        elif arg in inputs:
            return inputs[arg]
        elif isinstance(arg, InputSelector):
            return inputs[arg.input_to_select]
        elif isinstance(arg, InputNodeSelector):
            return inputs[arg.id()]

        if isinstance(arg, Operator):
            return arg.evaluate(inputs)

        raise Exception(f"Janis cannot evaluate '{arg.__class__.__name__}'")

    def rewrite_operator(self, args_to_rewrite: dict):
        return self.__class__(*self.substitute_arg(args_to_rewrite, self.args))

    @staticmethod
    def substitute_arg(args_to_rewrite: dict, arg: any):
        if isinstance(arg, list):
            return [Operator.substitute_arg(args_to_rewrite, a) for a in arg]
        elif isinstance(arg, dict):
            return {
                k: Operator.substitute_arg(args_to_rewrite, a) for k, a in arg.items()
            }
        if arg in args_to_rewrite:
            return args_to_rewrite[arg]
        return arg

    @abstractmethod
    def evaluate(self, inputs):
        # I need each operator to define this method
        pass

    @abstractmethod
    def to_wdl(self, unwrap_operator, *args):
        pass

    @abstractmethod
    def to_cwl(self, unwrap_operator, *args):
        pass

    def to_string_formatter(self):
        import re
        from janis_core.operators.stringformatter import StringFormatter

        key = re.sub(r"\W+", "", str(self))
        kwarg = {key: self}
        return StringFormatter(f"{{{key}}}", **kwarg)

    def __repr__(self):
        return f"{self.__class__.__name__}({', '.join(repr(a) for a in self.args)})"


class IndexOperator(Operator, ABC):
    def __init__(self, base, index):
        super().__init__(base, index)

    @staticmethod
    def friendly_signature():
        return "Array[X], Int -> X"

    def argtypes(self):
        return [Array(AnyType), Int]

    def returntype(self):
        return self.args[0].returntype().subtype()

    def __str__(self):
        base, index = self.args
        return f"{base}[{index}]"

    def __repr__(self):
        return str(self)

    def evaluate(self, inputs):
        iterable, idx = self.evaluate_arg(self.args, inputs)

        return iterable[idx]

    def to_wdl(self, unwrap_operator, *args):
        base, index = [unwrap_operator(a) for a in self.args]
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

    @staticmethod
    @abstractmethod
    def apply_to(value):
        pass

    def __str__(self):
        internal = self.args[0]
        return f"{self.symbol()}({internal})"

    def __repr__(self):
        return str(self)

    def evaluate(self, inputs):
        result = self.evaluate_arg(self.args[0], inputs)
        return self.apply_to(result)

    def to_wdl(self, unwrap_operator, *args):
        return f"{self.wdl_symbol()}({unwrap_operator(*args)})"

    def to_cwl(self, unwrap_operator, *args):
        return f"{self.cwl_symbol()}({unwrap_operator(*args)})"


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

    @staticmethod
    @abstractmethod
    def apply_to(arg1, arg2):
        pass

    def evaluate(self, inputs):
        arg1, arg2 = [self.evaluate_arg(a, inputs) for a in self.args]
        return self.apply_to(arg1, arg2)

    def to_wdl(self, unwrap_operator, *args):
        arg1, arg2 = [unwrap_operator(a) for a in self.args]
        return f"({arg1} {self.wdl_symbol()} {arg2})"

    def to_cwl(self, unwrap_operator, *args):
        arg1, arg2 = [unwrap_operator(a) for a in self.args]
        return f"({arg1} {self.cwl_symbol()} {arg2})"

    def __str__(self):
        args = self.args
        return f"({args[0]} {self.symbol()} {args[1]})"

    def __repr__(self):
        return str(self)


class AsStringOperator(SingleValueOperator):
    @staticmethod
    def symbol():
        return "str"

    @staticmethod
    def friendly_signature():
        return "X -> String"

    @staticmethod
    def wdl_symbol():
        return ""

    @staticmethod
    def cwl_symbol():
        return "String"

    @staticmethod
    def apply_to(value):
        return str(value)

    def argtypes(self):
        return [AnyType]

    def returntype(self):
        return String


class AsBoolOperator(SingleValueOperator):
    @staticmethod
    def symbol():
        return "bool"

    @staticmethod
    def friendly_signature():
        return "X -> Boolean"

    @staticmethod
    def wdl_symbol():
        return ""

    @staticmethod
    def cwl_symbol():
        return "Boolean"

    @staticmethod
    def apply_to(value):
        return bool(value)

    def argtypes(self):
        return [AnyType]

    def returntype(self):
        return Boolean


class AsIntOperator(SingleValueOperator):
    @staticmethod
    def symbol():
        return "int"

    @staticmethod
    def friendly_signature():
        return "X -> Int"

    @staticmethod
    def wdl_symbol():
        return ""

    @staticmethod
    def cwl_symbol():
        return "Number"

    def argtypes(self):
        return [AnyType]

    def returntype(self):
        return Int

    @staticmethod
    def apply_to(value):
        return int(value)


class AsFloatOperator(SingleValueOperator):
    @staticmethod
    def symbol():
        return "float"

    @staticmethod
    def friendly_signature():
        return "X -> Float"

    @staticmethod
    def wdl_symbol():
        return ""

    @staticmethod
    def cwl_symbol():
        return "Number"

    @staticmethod
    def apply_to(value):
        return float(value)

    def argtypes(self):
        return [AnyType]

    def returntype(self):
        return Float()
