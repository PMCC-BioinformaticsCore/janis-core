from copy import copy
from typing import List

from janis_core.utils.logger import Logger
from janis_core.types import (
    DataType,
    UnionType,
    Int,
    File,
    Float,
    Directory,
    get_instantiated_type,
)

from janis_core.types.common_data_types import String, Array, AnyType
from janis_core.operators.operator import Operator, InputSelector


class ReadContents(Operator):
    @staticmethod
    def friendly_signature():
        return "File -> String"

    def argtypes(self) -> List[DataType]:
        return [File()]

    def to_python(self, unwrap_operator, *args):
        raise NotImplementedError("Determine _safe_ one line solution for ReadContents")

    def to_wdl(self, unwrap_operator, *args):
        arg = unwrap_operator(args[0])
        return f"read_string({arg})"

    def to_cwl(self, unwrap_operator, *args):
        arg = unwrap_operator(
            args[0], add_path_suffix_to_input_selector_if_required=False
        )
        return f"{arg}.contents"

    def returntype(self):
        return String()

    def requires_contents(self):
        return True

    def evaluate(self, inputs):
        file = self.evaluate_arg(self.args[0], inputs)
        with open(file) as f:
            return f.read()


class ReadJsonOperator(Operator):
    @staticmethod
    def friendly_signature():
        return f"File -> Dict[str, any]"

    def evaluate(self, inputs):
        file = self.evaluate_arg(self.args[0], inputs)
        from json import load

        with open(file) as f:
            return load(f)

    def to_python(self, unwrap_operator, *args):
        raise NotImplementedError("Determine _safe_ one line solution for ReadContents")

    def to_wdl(self, unwrap_operator, *args):
        f = unwrap_operator(self.args[0])
        return f"read_json({f})"

    def to_cwl(self, unwrap_operator, *args):
        fp = unwrap_operator(
            self.args[0], add_path_suffix_to_input_selector_if_required=False
        )
        return f"JSON.parse({fp}.contents)"

    def requires_contents(self):
        return True

    def argtypes(self) -> List[DataType]:
        return [File()]

    def returntype(self):
        # dictionary?
        return String


class JoinOperator(Operator):
    def __init__(self, iterable, separator):
        super().__init__(iterable, separator)

    @staticmethod
    def friendly_signature():
        return "Array[X], String -> String"

    def argtypes(self):
        return [Array(UnionType(AnyType)), String]

    def returntype(self):
        return String()

    def to_python(self, unwrap_operator, *args):
        iterable, separator = [unwrap_operator(a) for a in self.args]
        return f"{separator}.join({iterable})"

    def to_wdl(self, unwrap_operator, *args):
        iterable, separator = [unwrap_operator(a) for a in self.args]
        iterable_arg = self.args[0]
        if isinstance(iterable_arg, list):
            is_optional = any(
                get_instantiated_type(a.returntype()).optional for a in iterable_arg
            )
        else:
            rettype = get_instantiated_type(iterable_arg.returntype())
            if rettype.is_array():
                is_optional = rettype.subtype().optional
            else:
                is_optional = rettype.optional

        if is_optional:
            return f"sep({separator}, select_first([{iterable}, []]))"
        else:
            return f"sep({separator}, {iterable})"

    def to_cwl(self, unwrap_operator, *args):
        iterable, separator = [unwrap_operator(a) for a in self.args]
        return f"{iterable}.join({separator})"

    def evaluate(self, inputs):
        iterable, separator = self.evaluate_arg(self.args, inputs)
        return str(separator).join((str(el) for el in iterable))


class BasenameOperator(Operator):
    @staticmethod
    def friendly_signature():
        return "Union[File, Directory] -> String"

    def to_python(self, unwrap_operator, *args):
        arg = unwrap_operator(args[0])
        return f"os.path.basename({arg})"

    def to_wdl(self, unwrap_operator, *args):
        arg = unwrap_operator(args[0])
        return f"basename({arg})"

    def to_cwl(self, unwrap_operator, *args):
        arg = unwrap_operator(
            args[0], add_path_suffix_to_input_selector_if_required=False
        )
        return arg + ".basename"

    def argtypes(self):
        return [UnionType(File, Directory)]

    def returntype(self):
        return String()

    def __str__(self):
        return str(self.args[0]) + ".basename"

    def __repr__(self):
        return str(self)

    def evaluate(self, inputs):
        from os.path import basename

        return basename(self.evaluate_arg(self.args[0], inputs))


class TransposeOperator(Operator):
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

    def to_python(self, unwrap_operator, *args):
        iterable = unwrap_operator(self.args[0])
        return f"[[{iterable}[j][i] for j in range(len({iterable}))] for i in range(len({iterable}[0]))]"

    def to_wdl(self, unwrap_operator, *args):
        return f"transform({unwrap_operator(args[0])})"

    def to_cwl(self, unwrap_operator, *args):
        return (
            unwrap_operator(args[0])
            + ".reduce(function(prev, next) { return next.map(function(item, i) { return (prev[i] || []).concat(next[i]); }) }, [])"
        )

    def evaluate(self, inputs):
        ar = self.evaluate_arg(self.args[0], inputs)
        return [[ar[i][j] for i in range(len(ar))] for j in range(len(ar[0]))]


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

    def to_python(self, unwrap_operator, *args):
        arg = unwrap_operator(self.args[0])
        return f"len({arg})"

    def to_wdl(self, unwrap_operator, *args):
        arg = unwrap_operator(self.args[0])
        return f"length({arg})"

    def to_cwl(self, unwrap_operator, *args):
        arg = unwrap_operator(self.args[0])
        return f"{arg}.length"

    def evaluate(self, inputs):
        ar = self.evaluate_arg(self.args[0], inputs)
        return len(ar)


class RangeOperator(Operator):
    @staticmethod
    def friendly_signature():
        return "Int -> Array[Int]"

    def argtypes(self):
        return [Int]

    def returntype(self):
        return Array(Int())

    def __str__(self):
        return f"0...{self.args[0]}"

    def __repr__(self):
        return str(self)

    def to_python(self, unwrap_operator, *args):
        arg = unwrap_operator(self.args[0])
        return f"range({arg})"

    def to_wdl(self, unwrap_operator, *args):
        arg = unwrap_operator(self.args[0])
        return f"range({arg})"

    def to_cwl(self, unwrap_operator, *args):
        arg = unwrap_operator(self.args[0])
        return f"Array.from({{ length: {arg} + 1 }}, (_, i) => i)"
        # return f"{arg}.length"

    def evaluate(self, inputs):
        ar = self.evaluate_arg(self.args[0], inputs)
        return list(range(ar))


class FlattenOperator(Operator):
    @staticmethod
    def friendly_signature():
        return "Array[Array[X]] -> Array[x]"

    def argtypes(self):
        return [Array(Array(AnyType))]

    def returntype(self):
        return Array(self.args[0].returntype().subtype().subtype())

    def __str__(self):
        return f"flatten({self.args[0]})"

    def __repr__(self):
        return str(self)

    def to_python(self, unwrap_operator, *args):
        arg = unwrap_operator(self.args[0])
        return f"[el for sublist in {arg} for el in sublist]"

    def to_wdl(self, unwrap_operator, *args):
        arg = unwrap_operator(self.args[0])
        return f"flatten({arg})"

    def to_cwl(self, unwrap_operator, *args):
        arg = unwrap_operator(self.args[0])
        return f"{arg}.flat()"

    def evaluate(self, inputs):
        ar = self.evaluate_arg(self.args[0], inputs)
        return [el for sl in ar for el in sl]


class ApplyPrefixOperator(Operator):
    @staticmethod
    def friendly_signature():
        return "String, Array[String] -> Array[String]"

    def argtypes(self):
        return [String, Array(AnyType)]

    def returntype(self):
        return Array(String)

    def __str__(self):
        prefix, iterable = self.args
        return f"['{prefix}' * {iterable}]"

    def __repr__(self):
        return str(self)

    def to_python(self, unwrap_operator, *args):
        prefix, iterable = [unwrap_operator(a) for a in self.args]
        return f"[{prefix} + i for i in {iterable}]"

    def to_wdl(self, unwrap_operator, *args):
        prefix, iterable = [unwrap_operator(a) for a in self.args]
        return f"prefix({prefix}, {iterable})"

    def to_cwl(self, unwrap_operator, *args):
        prefix, iterable = [unwrap_operator(a) for a in self.args]
        return f"{iterable}.map(function (inner) {{ return {prefix} + inner; }})"

    def evaluate(self, inputs):
        prefix, iterable = self.evaluate_arg(self.args, inputs)
        return [f"{prefix}{el}" for el in iterable]


class FileSizeOperator(Operator):
    """
    Returned in MB: Note that this does NOT include the reference files (yet)
    """

    def __new__(cls, *args, **kwargs):
        multiplier = None
        src, *otherargs = args

        if len(otherargs) == 1:
            f = otherargs[0].lower()
            multiplier_heirarchy = [
                ("ki" in f, 1024),
                ("k" in f, 1000),
                ("mi" in f, 1.024),
                ("gi" in f, 0.001024),
                ("g" in f, 0.001),
            ]
            if not any(m[0] for m in multiplier_heirarchy):
                Logger.warn(
                    f"Couldn't determine prefix {f} for FileSizeOperator, defaulting to MB"
                )
            else:
                multiplier = [m[1] for m in multiplier_heirarchy if m[0] is True][0]

        instance = super(FileSizeOperator, cls).__new__(cls)
        instance.__init__(args[0])

        if multiplier is not None and multiplier != 1:
            return instance * multiplier
        return instance

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

    def to_python(self, unwrap_operator, *args):
        f = unwrap_operator(self.args[0])
        return f"os.stat({f}).st_size / 1000"

    def to_wdl(self, unwrap_operator, *args):
        f = unwrap_operator(self.args[0])
        return f'size({f}, "MB")'

    def to_cwl(self, unwrap_operator, *args):
        f = unwrap_operator(
            self.args[0], add_path_suffix_to_input_selector_if_required=False
        )
        return f"({f}.size / 1048576)"

    def evaluate(self, inputs):
        from os.path import getsize

        file = self.evaluate_arg(self.args[0], inputs)
        return getsize(file) / 1048576


class FirstOperator(Operator):
    @staticmethod
    def friendly_signature():
        return "Array[X?] -> X"

    def argtypes(self):
        return [Array(AnyType)]

    def returntype(self):
        if isinstance(self.args[0], list):
            rettype = self.args[0][0].returntype()
        else:
            rettype = self.args[0].subtype()

        rettype = copy(rettype)
        rettype.optional = False
        return rettype

    def __str__(self):
        iterable = self.args[0]
        return f"first({iterable}]"

    def __repr__(self):
        return str(self)

    def to_python(self, unwrap_operator, *args):
        iterable = unwrap_operator(self.args[0])
        return f"[a for a in {iterable} if a is not None][0]"

    def to_wdl(self, unwrap_operator, *args):
        iterable = unwrap_operator(self.args[0])
        return f"select_first({iterable})"

    def to_cwl(self, unwrap_operator, *args):
        iterable = unwrap_operator(self.args[0])
        return f"{iterable}.filter(function (inner) {{ return inner != null }})[0]"

    def evaluate(self, inputs):
        iterable = self.evaluate_arg(self.args[0], inputs)
        return [i for i in iterable if i is not None][0]


class FilterNullOperator(Operator):
    @staticmethod
    def friendly_signature():
        return "Array[X?] -> Array[X]"

    def argtypes(self):
        return [Array(AnyType)]

    def returntype(self):
        if isinstance(self.args[0], list):
            rettype = self.args[0][0].returntype()
        else:
            outer_rettype = get_instantiated_type(self.args[0].returntype())
            if not isinstance(outer_rettype, Array):
                # hmmm, this could be a bad input selector
                rettype = outer_rettype
                if not isinstance(self.args[0], InputSelector):
                    Logger.warn(
                        f'Expected return type of "{self.args[0]}" to be an array, '
                        f'but found {outer_rettype}, will return this as a returntype.'
                    )
            else:
                rettype = outer_rettype.subtype()

        rettype = copy(get_instantiated_type(rettype))
        rettype.optional = False
        return Array(rettype)

    def __str__(self):
        iterable = self.args[0]
        return f"filter_null({iterable}]"

    def __repr__(self):
        return str(self)

    def to_python(self, unwrap_operator, *args):
        iterable = unwrap_operator(self.args[0])
        return f"[a for a in {iterable} if a is not None]"

    def to_wdl(self, unwrap_operator, *args):
        iterable = unwrap_operator(self.args[0])
        return f"select_all({iterable})"

    def to_cwl(self, unwrap_operator, *args):
        iterable = unwrap_operator(self.args[0])
        return f"{iterable}.filter(function (inner) {{ return inner != null }})"

    def evaluate(self, inputs):
        iterable = self.evaluate_arg(self.args[0], inputs)
        return [i for i in iterable if i is not None]


class ReplaceOperator(Operator):

    @staticmethod
    def friendly_signature():
        return "Base: String, Pattern: String, Replacement: String -> String"

    def argtypes(self) -> List[DataType]:
        return [String(), String(), String()]

    def evaluate(self, inputs):
        base, pattern, replacement = [self.evaluate_arg(a, inputs) for a in self.args]
        import re
        return re.sub(pattern, replacement, base)

    def to_wdl(self, unwrap_operator, *args):
        base, pattern, replacement = [unwrap_operator(a) for a in self.args]
        return f"sub({base}, {pattern}, {replacement})"

    def to_cwl(self, unwrap_operator, *args):
        base, pattern, replacement = [unwrap_operator(a) for a in self.args]
        return f"{base}.replace(new RegExp({pattern}), {replacement})"

    def to_python(self, unwrap_operator, *args):
        base, pattern, replacement = [unwrap_operator(a) for a in self.args]
        return f"re.sub({pattern}, {replacement}, {base})"

    def returntype(self) -> DataType:
        return String()


# class Stdout(Operator):
#     @staticmethod
#     def friendly_signature():
#         return "() -> File"
#
#     def argtypes(self) -> List[DataType]:
#         return []
#
#     def evaluate(self, inputs):
#         return ""
#
#     def to_wdl(self, unwrap_operator, *args):
#         return "stdout()"
#
#     def to_cwl(self, unwrap_operator, *args):
#         return "self[0]"
#
#     def __init__(self, return_type=File):
#         super().__init__()
#         self.return_type = get_instantiated_type(return_type)
#
#     def returntype(self):
#         return self.return_type
#
#     def to_string_formatter(self):
#         kwargs = {"stdout": self}
#         from janis_core.operators.stringformatter import StringFormatter
#
#         return StringFormatter("{stdout}", **kwargs)
