from typing import List

from janis_core import DataType
from janis_core.types import UnionType, Int, File, Float, Directory

from janis_core.types.common_data_types import String, Array, AnyType
from janis_core.operators.operator import Operator


class ReadContents(Operator):
    @staticmethod
    def friendly_signature():
        return "File -> String"

    def argtypes(self) -> List[DataType]:
        return [File()]

    def to_wdl(self, unwrap_operator, *args):
        arg = unwrap_operator(args[0])
        return f"read_string({arg})"

    def to_cwl(self, unwrap_operator, *args):
        arg = unwrap_operator(args[0])
        return f"{arg}.contents"

    def returntype(self):
        return String

    def requires_contents(self):
        return True


class JoinOperator(Operator):
    @staticmethod
    def friendly_signature():
        return "Array[X], String -> String"

    def argtypes(self):
        return [Array(UnionType(AnyType)), String]

    def returntype(self):
        return String

    def to_wdl(self, unwrap_operator, *args):
        iterable, separator = [unwrap_operator(a) for a in self.args]
        return f"sep({separator}, select_first([{iterable}, []]))"

    def to_cwl(self, unwrap_operator, *args):
        iterable, separator = [unwrap_operator(a) for a in self.args]
        return f"{iterable}.join({separator})"


class BasenameOperator(Operator):
    @staticmethod
    def friendly_signature():
        return "Union[File, Directory] -> String"

    def to_wdl(self, unwrap_operator, *args):
        arg = args[0]
        return f"basename({unwrap_operator(arg)})"

    def to_cwl(self, unwrap_operator, *args):
        return unwrap_operator(args[0]) + ".basename"

    def argtypes(self):
        return [UnionType(File, Directory)]

    def returntype(self):
        return String()

    def __str__(self):
        return str(self.args[0]) + ".basename"

    def __repr__(self):
        return str(self)


class TransformOperator(Operator):
    @staticmethod
    def friendly_signature():
        return "Array[Array[X]] -> Array[Array[x]]"

    def argtypes(self):
        return [Array(Array(AnyType))]

    def returntype(self):
        return Array(Array(AnyType))

    def __str__(self):
        return str(f"transform({self.args[0]})")

    def __repr__(self):
        return str(self)

    def to_wdl(self, unwrap_operator, *args):
        return f"transform({unwrap_operator(args[0])})"

    def to_cwl(self, unwrap_operator, *args):
        return (
            unwrap_operator(args[0])
            + ".map(function(c, i) { return q.map(function(row) { return row[i]; }); })"
        )


class LengthOperator(Operator):
    @staticmethod
    def friendly_signature():
        return "Array[X] -> Int"

    def argtypes(self):
        return [Array(AnyType)]

    def returntype(self):
        return Int()

    def __str__(self):
        return f"{self.args[0]}.length"

    def __repr__(self):
        return str(self)

    def to_wdl(self, unwrap_operator, *args):
        arg = unwrap_operator(self.args[0])
        return f"length({arg})"

    def to_cwl(self, unwrap_operator, *args):
        arg = unwrap_operator(self.args[0])
        return f"{arg}.length"


class FlattenOperator(Operator):
    @staticmethod
    def friendly_signature():
        return "Array[Array[X]] -> Array[x]"

    def argtypes(self):
        return [Array(Array(AnyType))]

    def returntype(self):
        return Array(self.args[0].subtype().subtype())

    def __str__(self):
        return f"flatten({self.args[0]})"

    def __repr__(self):
        return str(self)

    def to_wdl(self, unwrap_operator, *args):
        arg = unwrap_operator(self.args[0])
        return f"flatten({arg})"

    def to_cwl(self, unwrap_operator, *args):
        arg = unwrap_operator(self.args[0])
        return f"{arg}.flat()"


class ApplyPrefixOperator(Operator):
    @staticmethod
    def friendly_signature():
        return "String, Array[String] -> String"

    def argtypes(self):
        return [String, Array(AnyType)]

    def returntype(self):
        return Array(String)

    def __str__(self):
        prefix, iterable = self.args
        return f"['{prefix}' * {iterable}]"

    def __repr__(self):
        return str(self)

    def to_wdl(self, unwrap_operator, *args):
        prefix, iterable = [unwrap_operator(a) for a in self.args]
        return f"prefix({prefix}, {iterable})"

    def to_cwl(self, unwrap_operator, *args):
        prefix, iterable = [unwrap_operator(a) for a in self.args]
        return f"{iterable}.map(function (inner) {{ return {prefix} + inner; }})"


class FileSizeOperator(Operator):
    """
    Returned in MB: Note that this does NOT include the reference files (yet)
    """

    @staticmethod
    def friendly_signature():
        return "File -> Float"

    def argtypes(self):
        return [File()]

    def returntype(self):
        return Float

    def __str__(self):
        f = self.args[0]
        return f"file_size({f})"

    def __repr__(self):
        return str(self)

    def to_wdl(self, unwrap_operator, *args):
        f = unwrap_operator(self.args[0])
        return f'size({f}, "MB")'

    def to_cwl(self, unwrap_operator, *args):
        f = unwrap_operator(self.args[0])
        return f"({f}.size / 1048576)"


class FirstOperator(Operator):
    @staticmethod
    def friendly_signature():
        return "Array[X?] -> X"

    def argtypes(self):
        return [Array(AnyType)]

    def returntype(self):
        return self.args[0].subtype()

    def __str__(self):
        iterable = self.args[0]
        return f"first({iterable}]"

    def __repr__(self):
        return str(self)

    def to_wdl(self, unwrap_operator, *args):
        iterable = unwrap_operator(self.args[0])
        return f"select_first({iterable})"

    def to_cwl(self, unwrap_operator, *args):
        iterable = unwrap_operator(self.args[0])
        return f"{iterable}.filter(function (inner) {{ return inner != null }})[0]"


class FilterNullOperator(Operator):
    @staticmethod
    def friendly_signature():
        return "Array[X?] -> Array[X]"

    def argtypes(self):
        return [Array(AnyType)]

    def returntype(self):
        return self.args[0].subtype()

    def __str__(self):
        iterable = self.args[0]
        return f"filter_null({iterable}]"

    def __repr__(self):
        return str(self)

    def to_wdl(self, unwrap_operator, *args):
        iterable = unwrap_operator(self.args[0])
        return f"select_all({iterable})"

    def to_cwl(self, unwrap_operator, *args):
        iterable = unwrap_operator(self.args[0])
        return f"{iterable}.filter(function (inner) {{ return inner != null }})"
