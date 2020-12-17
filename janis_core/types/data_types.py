"""
    Custom Types::

    We want to define some standard interface that allows for the implementation of types.
    Types must be supported as Inputs and Outputs, and must be comparable, ie we can compare
    the output of one step to another step. We'll guarantee that

    We'll also may require separate interfaces for CWL and WDL implementations, we may not be
    able to genericise the interface enough to support an automatic conversion

    We are allowed to require that a type must register itself when loaded.

"""
from abc import ABC, abstractmethod
from typing import Any, List, Optional, Union, Type

from janis_core.deps import cwlgen, wdlgen

from janis_core.utils import is_array_prefix
from janis_core.utils.logger import Logger

NativeType = str
PythonPrimitive = Union[str, float, int, bool]


def is_python_primitive(t):
    return (
        isinstance(t, list) and len(t) > 0 and is_python_primitive(t[0])
    ) or isinstance(t, (str, float, int, bool))


# see below for ParseableType


class NativeTypes:
    kStr: NativeType = "str"
    kInt: NativeType = "int"
    kLong: NativeType = "long"
    kFloat: NativeType = "float"
    kBool: NativeType = "bool"
    kDouble: NativeType = "double"
    kFile: NativeType = "file"
    kDirectory: NativeType = "dir"
    kArray: NativeType = "array"
    kStdout: NativeType = "stdout"
    kStderr: NativeType = "stderr"

    _primitives: List[NativeType] = [kStr, kInt, kFloat, kLong, kDouble, kBool, kDouble]
    all: List[NativeType] = _primitives + [kFile, kDirectory, kArray, kStdout, kStderr]

    @staticmethod
    def is_primitive(t: NativeType) -> bool:
        return t in NativeTypes._primitives

    @staticmethod
    def is_valid(t: str):
        return t in NativeTypes.all

    @staticmethod
    def map_to_primitive(t: NativeType):
        if t == NativeTypes.kBool:
            return str
        elif t == NativeTypes.kInt:
            return int
        elif t == NativeTypes.kLong:
            return int
        elif t == NativeTypes.kFloat:
            return float
        elif t == NativeTypes.kDouble:
            return float
        elif t == NativeTypes.kStr:
            return str
        elif t == NativeTypes.kFile:
            return str
        elif t == NativeTypes.kDirectory:
            return str
        elif t == NativeTypes.kArray:
            return list
        elif t == NativeTypes.kStdout:
            return str
        elif t == NativeTypes.kStderr:
            return str
        raise Exception(
            f"Unhandled primitive type {t}, expected one of {', '.join(NativeTypes.all)}"
        )

    @staticmethod
    def map_to_cwl(t: NativeType):
        if t == NativeTypes.kBool:
            return "boolean"
        elif t == NativeTypes.kInt:
            return "int"
        elif t == NativeTypes.kLong:
            return "long"
        elif t == NativeTypes.kFloat:
            return "float"
        elif t == NativeTypes.kDouble:
            return "double"
        elif t == NativeTypes.kStr:
            return "string"
        elif t == NativeTypes.kFile:
            return "File"
        elif t == NativeTypes.kDirectory:
            return "Directory"
        elif t == NativeTypes.kArray:
            return "array"
        elif t == NativeTypes.kStdout:
            return "stdout"
        elif t == NativeTypes.kStderr:
            return "stderr"
        raise Exception(
            f"Unhandled primitive type {t}, expected one of {', '.join(NativeTypes.all)}"
        )

    @staticmethod
    def map_to_wdl(t: NativeType):
        from janis_core.deps import wdlgen as wdl

        if t == NativeTypes.kBool:
            return wdl.PrimitiveType.kBoolean
        elif t == NativeTypes.kInt:
            return wdl.PrimitiveType.kInt

        elif (
            t == NativeTypes.kLong
            or t == NativeTypes.kFloat
            or t == NativeTypes.kDouble
        ):
            return wdl.PrimitiveType.kFloat
        elif t == NativeTypes.kStr:
            return wdl.PrimitiveType.kString
        elif t == NativeTypes.kFile:
            return wdl.PrimitiveType.kFile
        elif t == NativeTypes.kStdout:
            return wdl.PrimitiveType.kFile
        elif t == NativeTypes.kStderr:
            return wdl.PrimitiveType.kFile
        elif t == NativeTypes.kDirectory:
            Logger.log(
                "Using data_type 'Directory' for wdl, this requires cromwell>=37 and language=development"
            )
            return wdl.PrimitiveType.kDirectory
        elif t == NativeTypes.kArray:
            return wdl.ArrayType.kArray
        raise Exception(
            f"Unhandled primitive type {t}, expected one of {', '.join(NativeTypes.all)}"
        )


class DataType(ABC):
    def __init__(self, optional=False):
        self.optional = optional if optional is not None else False
        self.is_prim = NativeTypes.is_primitive(self.primitive())

    def is_array(self):
        return False

    def __repr__(self):
        return self.id()

    @staticmethod
    @abstractmethod
    def name() -> str:
        raise Exception("Subclass MUST override name field")

    @classmethod
    def __hash__(cls):
        return cls.name()

    @staticmethod
    def secondary_files() -> Optional[List[str]]:
        return None

    @staticmethod
    @abstractmethod
    def primitive() -> NativeType:
        raise Exception("Subclass MUST override the 'primitive' method")

    @staticmethod
    @abstractmethod
    def doc() -> str:
        """
        Subclasses should override this class to provide additional information on how to
        correctly provide data to the class, what inputs it may have and what other types
        are compatible
        """
        raise Exception("Subclass MUST override the 'doc' field")

    # The following methods don't need to be overriden, but can be

    def id(self):
        if self.optional:
            return f"Optional<{self.name()}>"
        return self.name()

    @abstractmethod
    def validate_value(self, meta: Any, allow_null_if_not_optional: bool) -> bool:
        pass

    def coerce_value_if_possible(self, value):
        return value

    @abstractmethod
    def invalid_value_hint(self, meta):
        pass

    def identify(self):
        print(self.id())

    def can_receive_from(self, other, source_has_default=False) -> bool:
        """
        Can this class receive from $other, likely going to be type(a) == type(b)
        :param other:
        :param source_has_default: If the source has default, then we can return true even if the source is optional
        :return:
        """

        # Depending on the way types are imported, the 'isinstance' method doesn't always work
        #
        # A small example:
        #     from janis import File as F1
        #     from janis_core import File as F2
        #
        # Although these are the same definition, they won't actually compare to the same value

        if other.name().lower() == "union":
            return all(
                self.can_receive_from(t, source_has_default=source_has_default)
                for t in other.subtypes
            )

        receive_from = list(
            reversed([x.__name__ for x in type(other.received_type()).mro()])
        )
        receive_to = list(
            reversed([x.__name__ for x in type(self.received_type()).mro()])
        )

        if not is_array_prefix(receive_to, receive_from):
            return False
        if self.optional or source_has_default:
            # If I'm optional I can receive from optional / non optional
            return True
        # If I'm not optional, I must receive from not optional
        return not other.optional

    def received_type(self):
        """
        The type that will be received if joined from this type
        :return: mostly self, except for STDOUT | STDERR | STDIN
        """
        return self

    def input_field_from_input(self, meta):
        """
        Method to convert the field definition into a generic CWL-esque response
        :param meta:
        :return:
        """
        return None

    def _question_mark_if_optional(self, has_default: bool = False):
        return "?" if self.optional or has_default else ""

    def cwl_type(
        self, has_default=False
    ) -> Union[str, cwlgen.Type, List[Union[str, cwlgen.Type]]]:
        tp = NativeTypes.map_to_cwl(self.primitive())
        return (
            [tp, "null"] if self.optional and not has_default else tp
        )  # and not has_default

    def map_cwl_type(self, parameter: cwlgen.Parameter) -> cwlgen.Parameter:
        if not NativeTypes.is_valid(self.primitive()):
            raise Exception(
                f"{self.id()} must declare its primitive as one of the NativeTypes "
                f"({', '.join(NativeTypes.all)})"
            )

        tp = NativeTypes.map_to_cwl(self.primitive())
        parameter.type = [tp, "null"] if self.optional else tp
        parameter.secondaryFiles = self.secondary_files()
        return parameter

    def cwl_input(self, value: Any):
        return value

    def wdl(self, has_default=False) -> wdlgen.WdlType:
        qm = self._question_mark_if_optional(has_default)
        return wdlgen.WdlType.parse_type(NativeTypes.map_to_wdl(self.primitive()) + qm)

    def parse_value(self, valuetoparse):
        """
        Sometimes it's useful for a value to be parsed if possible by this class.
        Eg, number, or arrays.
        """
        return valuetoparse

    def copy(self):
        from copy import deepcopy

        return deepcopy(self)

    # def default(self):
    #     return self.default_value


ParseableTypeBase = Union[Type[PythonPrimitive], DataType, Type[DataType]]
ParseableType = ParseableTypeBase
