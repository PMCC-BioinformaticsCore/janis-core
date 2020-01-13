###################
# Implementations #
###################
from inspect import isclass
from typing import Union, Type, Dict, Any

import cwlgen
import wdlgen

from janis_core.types.data_types import (
    DataType,
    NativeTypes,
    NativeType,
    PythonPrimitive,
)
from janis_core.utils.generics_util import is_generic, is_qualified_generic

ParseableTypeBase = Union[Type[PythonPrimitive], DataType, Type[DataType]]
ParseableType = ParseableTypeBase


class String(DataType):
    @staticmethod
    def name():
        return "String"

    @staticmethod
    def primitive():
        return NativeTypes.kStr

    @staticmethod
    def doc():
        return "A string"

    @classmethod
    def schema(cls) -> Dict:
        return {"type": "string", "required": True}

    def input_field_from_input(self, meta):
        return next(iter(meta.values()))

    def can_receive_from(self, other, source_has_default=False):
        if isinstance(other, Filename):
            return True
        return super().can_receive_from(other, source_has_default=source_has_default)

    def validate_value(self, meta: Any, allow_null_if_not_optional: bool):
        if meta is None:
            return self.optional or allow_null_if_not_optional
        return isinstance(meta, str)

    def invalid_value_hint(self, meta):
        if meta is None:
            return "value was null"
        if self.validate_value(meta, True):
            return None

        return f"Value was of type {type(meta)}, expected string"


class Filename(String):
    def __init__(
        self, prefix="generated", suffix=None, extension: str = None, guid: str = None
    ):
        """
        :param suffix: suffix the guid
        :param extension: with no '.' (dot)
        :param guid: Use this guid instead of generating one
        """
        import uuid

        self.prefix = prefix
        self.extension = extension
        self.suffix = suffix
        self.guid = guid if guid is not None else str(uuid.uuid1())
        super().__init__(optional=True)

    @staticmethod
    def name() -> str:
        return "Filename"

    @staticmethod
    def primitive() -> NativeType:
        return NativeTypes.kStr

    def cwl_type(self, has_default=False):
        self.optional = False
        t = super().cwl_type()
        self.optional = True
        return t

    @staticmethod
    def doc() -> str:
        return """
This class is a placeholder for generated filenames, by default it is optional and CAN be overrided, 
however the program has been structured in a way such that these names will be generated based on the step label. 
These should only be used when the tool _requires_ a filename to output and you aren't 
concerned what the filename should be. The Filename DataType should NOT be used as an output.
""".strip()

    @classmethod
    def schema(cls) -> Dict:
        pass

    def map_cwl_type(self, parameter: cwlgen.Parameter):
        super().map_cwl_type(parameter)
        parameter.default = self.generated_filenamecwl()

    def generated_filename(self) -> str:

        pre = self.prefix
        suf = ("-" + str(self.suffix)) if self.suffix else ""
        ex = "" if self.extension is None else self.extension
        return pre + suf + ex

    def generated_filenamecwl(self) -> str:
        return f'"{self.generated_filename()}"'
        # code = "Math.random().toString(16).substring(2, 8)"
        # pf = (self.prefix + "-") if self.prefix else ""
        # sf = self.suffix if self.suffix else ""
        # ext = self.extension if self.extension else ""
        # return f'"{pf}generated-" + {code} + "{sf + ext}"'

    def can_receive_from(self, other: DataType, source_has_default=False):
        # Specific override because Filename should be able to receive from string
        if isinstance(other, String):
            return True  # Always provides default, and is always optional
        return super().can_receive_from(other, source_has_default=source_has_default)

    def wdl(self, has_default=True):
        return wdlgen.WdlType.parse_type(NativeTypes.map_to_wdl(self.primitive()))

    def validate_value(self, meta: Any, allow_null_if_not_optional: bool):
        return True

    def invalid_value_hint(self, meta):
        return None


class Int(DataType):
    @staticmethod
    def name():
        return "Integer"

    @staticmethod
    def primitive():
        return NativeTypes.kInt

    def doc(self):
        return "An integer"

    @classmethod
    def schema(cls) -> Dict:
        return {"type": "number", "required": True}

    def input_field_from_input(self, meta):
        return next(iter(meta.values()))

    def validate_value(self, meta: Any, allow_null_if_not_optional: bool):
        if meta is None:
            return self.optional or allow_null_if_not_optional
        return isinstance(meta, int)

    def invalid_value_hint(self, meta):
        if meta is None:
            return "value was null"
        if self.validate_value(meta, True):
            return None
        return f"Value was of type {type(meta)}, expected int"


class Float(DataType):
    @staticmethod
    def name():
        return "Float"

    @staticmethod
    def primitive():
        return NativeTypes.kFloat

    def doc(self):
        return "A float"

    @classmethod
    def schema(cls) -> Dict:
        return {"type": "number", "required": True}

    def input_field_from_input(self, meta: Dict):
        return next(iter(meta.values()))

    def validate_value(self, meta: Any, allow_null_if_not_optional: bool) -> bool:
        if meta is None:
            return self.optional or allow_null_if_not_optional
        return isinstance(meta, float) or isinstance(meta, int)

    def invalid_value_hint(self, meta):
        if meta is None:
            return "value was null"
        if self.validate_value(meta, True):
            return None
        return f"Value was of type {type(meta)}, expected float | int"


class Double(DataType):
    @staticmethod
    def name():
        return "Double"

    @staticmethod
    def primitive():
        return NativeTypes.kDouble

    def doc(self):
        return "An integer"

    @classmethod
    def schema(cls) -> Dict:
        return {"type": "number", "required": True}

    def input_field_from_input(self, meta):
        return next(iter(meta.values()))

    def validate_value(self, meta: Any, allow_null_if_not_optional: bool) -> bool:
        if meta is None:
            return self.optional or allow_null_if_not_optional
        return isinstance(meta, float) or isinstance(meta, int)

    def invalid_value_hint(self, meta):
        if meta is None:
            return "value was null"
        if self.validate_value(meta, True):
            return None
        return f"Value was of type {type(meta)}, expected float | int"


class Boolean(DataType):
    @staticmethod
    def name():
        return "Boolean"

    @staticmethod
    def primitive():
        return NativeTypes.kBool

    def doc(self):
        return "A boolean"

    @classmethod
    def schema(cls) -> Dict:
        return {"type": "boolean", "required": True}

    def input_field_from_input(self, meta):
        return next(iter(meta.values()))

    def validate_value(self, meta: Any, allow_null_if_not_optional: bool) -> bool:
        if meta is None:
            return self.optional or allow_null_if_not_optional
        return isinstance(meta, bool)

    def invalid_value_hint(self, meta):
        if meta is None:
            return "value was null"
        if self.validate_value(meta, True):
            return None
        return f"Value was of type {type(meta)}, expected bool"


class File(DataType):
    def __init__(self, optional=False, extension=None):
        """
        :param optional:
        :param common_extension: Used in CWL to try and guess the file extension where it's not available otherwise
        """
        super(File, self).__init__(optional=optional)
        self.extension = extension

    @staticmethod
    def name():
        return "File"

    @staticmethod
    def primitive():
        return NativeTypes.kFile

    def doc(self):
        return "A local file"

    @classmethod
    def schema(cls) -> Dict:
        return {"path": {"type": "string", "required": True}}

    def get_value_from_meta(self, meta):
        return meta.get("path")

    def cwl_input(self, value: Any):
        return {"class": cwlgen.CwlTypes.FILE, "path": value}

    def validate_value(self, meta: Any, allow_null_if_not_optional: bool) -> bool:
        if meta is None:
            return self.optional or allow_null_if_not_optional
        return isinstance(meta, str)

    def invalid_value_hint(self, meta):
        if meta is None:
            return "value was null"
        if self.validate_value(meta, True):
            return None

        return f"Value was of type {type(meta)}, expected string (path)"


class Directory(DataType):
    def __init__(self, optional=False):
        """
        Specifically exclude default
        """
        super(Directory, self).__init__(optional=optional)

    @staticmethod
    def name():
        return "Directory"

    @staticmethod
    def primitive():
        return NativeTypes.kDirectory

    def doc(self):
        return "A directory of files"

    def get_value_from_meta(self, meta):
        return meta["path"]

    @classmethod
    def schema(cls) -> Dict:
        return {"path": {"type": "string", "required": True}}

    def input_field_from_input(self, meta):
        return meta["path"]

    def cwl_input(self, value: Any):
        # WDL: "{workflowName}.label" = meta["path"}
        return {"class": cwlgen.CwlTypes.DIRECTORY, "path": value}

    def validate_value(self, meta: Any, allow_null_if_not_optional: bool) -> bool:
        if meta is None:
            return self.optional or allow_null_if_not_optional
        return isinstance(meta, str)

    def invalid_value_hint(self, meta):
        if meta is None:
            return "value was null"
        if self.validate_value(meta, True):
            return None
        return f"Value was of type {type(meta)}, expected string (path)"


class Array(DataType):
    def __init__(self, t: ParseableType, optional=False):
        resolvedtype = get_instantiated_type(t)
        if not isinstance(resolvedtype, DataType):
            raise Exception(f"Type t ({type(t)}) must be an instance of 'DataType'")

        self._t = resolvedtype
        super().__init__(optional)

    def subtype(self):
        return self._t

    @staticmethod
    def name():
        return "Array"

    @staticmethod
    def primitive():
        return NativeTypes.kArray

    def id(self):
        if self._t is None:
            return super().id()
        t = self._t
        typed = f"Array<{t.id()}>"
        if self.optional:
            return f"Optional<{typed}>"
        return typed

    def doc(self):
        return "An array"

    @classmethod
    def schema(cls) -> Dict:
        return {"type": "array"}

    def cwl_type(self, has_default=False):
        inp = cwlgen.CommandInputArraySchema(
            items=self._t.cwl_type(),
            # label=None,
            # input_binding=None
        )
        return [inp, "null"] if self.optional and not has_default else inp

    def map_cwl_type(self, parameter: cwlgen.Parameter) -> cwlgen.Parameter:
        parameter.type = cwlgen.CommandInputArraySchema(
            items=None, label=None, input_binding=None
        )
        return parameter

    def cwl_input(self, value: Any):
        if isinstance(value, list):
            return [self._t.cwl_input(v) for v in value]
        if value is None:
            return None
        else:
            raise Exception(f"Input value for input '{self.id()}' was not an array")

    def wdl(self, has_default=False) -> wdlgen.WdlType:
        ar = wdlgen.ArrayType(self._t.wdl(has_default=False), requires_multiple=False)
        return wdlgen.WdlType(ar, optional=self.optional or has_default)

    def can_receive_from(self, other, source_has_default=False):
        if isinstance(other, Array):
            return self._t.can_receive_from(other._t)
        if not self._t.can_receive_from(other):
            return False
        return super().can_receive_from(other, source_has_default=source_has_default)

    def input_field_from_input(self, meta):
        return next(iter(meta.values()))

    def validate_value(self, meta: Any, allow_null_if_not_optional: bool) -> bool:
        if meta is None:
            return self.optional or allow_null_if_not_optional
        if not isinstance(meta, list):
            return False
        return all(
            self.subtype().validate_value(q, allow_null_if_not_optional) for q in meta
        )

    def invalid_value_hint(self, meta):
        if meta is None:
            return "value was null"
        if self.validate_value(meta, True):
            return None

        if not isinstance(meta, list):
            return f"Value was of type {type(meta)}, expected type Array<{self.subtype().id()}>"

        hints = []
        st = self.subtype()
        for i in range(len(meta)):
            hint = st.invalid_value_hint(meta[i])
            if not hint:
                continue
            hints.append(f"{i}. {hint}")

        return str(hints)

    def parse_value(self, valuetoparse):
        if not isinstance(valuetoparse, list):
            valuetoparse = [valuetoparse]

        return [self.subtype().parse_value(v) for v in valuetoparse]

    def fundamental_type(self) -> DataType:
        st = self.subtype()
        if isinstance(st, Array):
            return st.fundamental_type()
        return st.received_type()


class Stdout(File):
    @staticmethod
    def name():
        return "Stdout"

    def __init__(self, subtype=None, stdoutname=None, optional=None):
        super().__init__(optional=False)

        subtype = get_instantiated_type(subtype) if subtype is not None else File()

        if subtype and not isinstance(subtype, File):
            raise Exception(
                "Janis does not currently support non-File stdout annotations"
            )

        self.stdoutname = stdoutname
        self.subtype = subtype

        if self.subtype.secondary_files():
            raise Exception(
                f"The subtype '{self.subtype.__name__}' has secondary files, "
                f"but stdout does not have the ability to collect files"
            )

    @staticmethod
    def primitive():
        return NativeTypes.kStdout

    def id(self):
        return f"stdout<{self.subtype.id()}>"

    def received_type(self):
        return self.subtype

    def validate_value(self, meta: Any, allow_null_if_not_optional: bool) -> bool:
        """
        Will always toss away the value
        """
        return True

    def invalid_value_hint(self, meta):
        return None


class Stderr(File):
    @staticmethod
    def name():
        return "Stderr"

    def __init__(self, subtype=None, stderrname=None):
        super().__init__(optional=False)

        subtype = get_instantiated_type(subtype) if subtype is not None else File()

        if subtype and not isinstance(subtype, File):
            raise Exception(
                "Janis does not currently support non-File stderr annotations"
            )

        self.stderrname = stderrname
        self.subtype = subtype

        if self.subtype.secondary_files():
            raise Exception(
                f"The subtype '{self.subtype.__name__}' has secondary files, "
                f"but stderr does not have the ability to collect files"
            )

    @staticmethod
    def primitive():
        return NativeTypes.kStderr

    def id(self):
        return f"stderr<{self.subtype.id()}>"

    def received_type(self):
        return self.subtype

    def validate_value(self, meta: Any, allow_null_if_not_optional: bool) -> bool:
        """
        Will always toss away the value
        """
        return True

    def invalid_value_hint(self, meta):
        return None


all_types = [
    String,
    Filename,
    Int,
    Float,
    Double,
    Boolean,
    File,
    Directory,
    Stdout,
    Array,
]


def get_from_python_type(dt, optional: bool = None, overrider=None):
    if dt is None:
        return None

    bc = overrider or get_instantiated_type
    typedt = type(dt)

    if dt == str or typedt == str:
        return String(optional=optional)
    if dt == bool or typedt == bool:
        return Boolean(optional=optional)
    if dt == int or typedt == int:
        return Int(optional=optional)
    if dt == float or typedt == float:
        return Float(optional=optional)

    if is_qualified_generic(dt):

        if str(dt).startswith("typing.List"):
            nt = bc(dt.__args__[0])
            return Array(nt, optional=optional)

        args = dt.__args__
        if len(args) > 2:
            raise Exception(f"Janis is unsure how to parse qualfied generic '{dt}'")

        aridxofnonetype = [
            i for i, val in enumerate(a == type(None) for a in args) if val
        ]
        optional = len(aridxofnonetype) > 0

        if len(aridxofnonetype) > 1 and optional is False:
            raise Exception("Janis cannot accept union ")

        idxofsubtype = (len(args) - 1 - aridxofnonetype[0]) if optional else 0
        subtype = args[idxofsubtype]

        nt = bc(subtype, optional=optional)
        return nt

    elif is_generic(dt):
        raise Exception(f"Generic {dt} was generic typing, but unqualified")


def get_instantiated_type(datatype: ParseableType, optional=None, overrider=None):

    bc = overrider or get_instantiated_type

    if isinstance(datatype, list):
        if len(datatype) == 0:
            raise TypeError("Couldn't determine type of array with length 0")
        return Array(bc(datatype[0]))

    if isinstance(datatype, DataType):
        return datatype

    if isclass(datatype) and issubclass(datatype, DataType):
        return datatype(optional=optional)

    dt = get_from_python_type(datatype, optional=optional, overrider=bc)
    if dt:
        return dt

    raise TypeError(f"Unable to parse type '{str(datatype)}'")
