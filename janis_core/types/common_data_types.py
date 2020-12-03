###################
# Implementations #
###################
from inspect import isclass
from typing import Dict, Any, Set


from janis_core.deps import cwlgen, wdlgen

from janis_core.utils.logger import Logger
from janis_core.__meta__ import GITHUB_URL
from janis_core.types.data_types import DataType, NativeTypes, NativeType, ParseableType
from janis_core.utils.generics_util import is_generic, is_qualified_generic


class UnionType(DataType):
    def __init__(self, *subtypes: ParseableType, optional=False):
        self._initial_subtypes = [t for t in subtypes]

        invalid_types = []
        valid_types = []
        types_with_secondaries = []
        for subtype in subtypes:
            resolvedtype = get_instantiated_type(subtype)
            if not isinstance(resolvedtype, DataType):
                invalid_types.append(resolvedtype)
            elif isinstance(resolvedtype, File) and resolvedtype.secondary_files():
                types_with_secondaries.append(types_with_secondaries)
            else:
                valid_types.append(resolvedtype)

        if len(types_with_secondaries) > 0:
            raise Exception(
                "UnionType doesn't accept data types with secondary files (yet), affected types: "
                + ", ".join(str(t) for t in types_with_secondaries)
            )

        if len(invalid_types) > 0:
            raise Exception(
                "UnionType contained invalid types "
                + ", ".join(str(t) for t in invalid_types)
            )

        if len(valid_types) < 1:
            raise Exception("UnionType is expecting at least 2 data type arguments")

        self.subtypes = valid_types
        super().__init__(optional)

    def is_array(self):
        return all(s.is_array() for s in self.subtypes)

    def id(self):
        return "Union<" + ", ".join(s.id() for s in self.subtypes) + ">"

    @staticmethod
    def name() -> str:
        return "Union"

    @staticmethod
    def primitive() -> NativeType:
        return None

    @staticmethod
    def doc() -> str:
        return "Union datatype"

    def validate_value(self, *args, **kwargs) -> bool:
        return any(t.validate_value(*args, **kwargs) for t in self.subtypes)

    def invalid_value_hint(self, *args, **kwargs):
        hints = [t.invalid_value_hint(*args, **kwargs) for t in self.subtypes]
        return ", ".join(t for t in hints if t)

    def can_receive_from(self, other, *args, **kwargs):
        if isinstance(other, UnionType):
            # we'll require all elements in the source to be received by this type-
            return all(
                self.can_receive_from(t, *args, **kwargs) for t in other.subtypes
            )
        return any(t.can_receive_from(other, *args, **kwargs) for t in self.subtypes)

    def wdl(self, has_default=False) -> wdlgen.WdlType:
        # custom stuff here
        wdl_data_types = [a.wdl() for a in self.subtypes]
        # we require the WDL to be identical for WDL to work

        if len(set(a.get_string() for a in wdl_data_types)) > 1:

            resuting_signatures = ", ".join(
                f"{a.id()}: {a.wdl().get_string()}" for a in self.subtypes
            )

            raise Exception(
                "Janis doesn't support UnionTypes in WDL where there is more than 1 WDL type signatures. "
                f"Please raise an issue on GitHub ({GITHUB_URL}) if this is a blocker. Resulting signatures: "
                + resuting_signatures
            )

        return wdl_data_types[0]

    def cwl_type(self, has_default=False):
        inner_types = [a.cwl_type(has_default=has_default) for a in self.subtypes]
        try:
            inner_types = list(set(inner_types))
        except Exception as e:
            Logger.debug(f"Error creating set from ({inner_types}): {e}")

        if len(inner_types) == 1:
            return inner_types[0]
        return inner_types


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
        return isinstance(meta, (str, float, int))

    def coerce_value_if_possible(self, value):
        return str(value)

    def invalid_value_hint(self, meta):
        if meta is None:
            return "value was null"
        if self.validate_value(meta, True):
            return None

        return f"Value was of type {type(meta)}, expected string"


class Filename(String):
    def __init__(
        self, prefix="generated", suffix=None, extension: str = None, optional=None
    ):
        """
        :param suffix: suffix the guid
        :param extension: with no '.' (dot)
        :param guid: Use this guid instead of generating one
        :param optional: IGNORED (legacy)
        """

        self.prefix = prefix
        self.extension = extension
        self.suffix = suffix

        super().__init__(optional=True)

    @staticmethod
    def name() -> str:
        return "Filename"

    @staticmethod
    def primitive() -> NativeType:
        return NativeTypes.kStr

    def cwl_type(self, has_default=False):
        return super().cwl_type()

    @staticmethod
    def doc() -> str:
        return """
This class is a placeholder for generated filenames, by default it is optional and CAN be overridden, 
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

    def generated_filename(self, replacements: Dict = None) -> str:
        repl = replacements or {}
        prefix = repl.get("prefix", self.prefix)
        suffix = repl.get("suffix", self.suffix)

        suf = ""
        if suffix:
            if str(suffix).startswith("."):
                suf = str(suffix)
            else:
                suf = "-" + str(suffix)
        ex = "" if self.extension is None else self.extension
        return prefix + suf + ex

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
        return super().wdl(has_default=has_default)

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

    def coerce_value_if_possible(self, value):
        try:
            return int(value)
        except:
            raise Exception(f"Value '{value}' cannot be coerced to an integer")

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

    def input_field_from_input(self, meta: Dict):
        return next(iter(meta.values()))

    def validate_value(self, meta: Any, allow_null_if_not_optional: bool) -> bool:
        if meta is None:
            return self.optional or allow_null_if_not_optional
        return isinstance(meta, float) or isinstance(meta, int)

    def coerce_value_if_possible(self, value):
        try:
            return float(value)
        except:
            raise Exception(f"Value '{value}' cannot be coerced to a float")

    def invalid_value_hint(self, meta):
        if meta is None:
            return "value was null"
        if self.validate_value(meta, True):
            return None
        return f"Value was of type {type(meta)}, expected float | int"


class Double(Float):
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

    def can_receive_from(self, other, *args, **kwargs) -> bool:
        if not other.optional and isinstance(other, Float):
            return True
        return super().can_receive_from(other, *args, **kwargs)


class Boolean(DataType):
    @staticmethod
    def name():
        return "Boolean"

    @staticmethod
    def primitive():
        return NativeTypes.kBool

    def doc(self):
        return "A boolean"

    def input_field_from_input(self, meta):
        return next(iter(meta.values()))

    def validate_value(self, meta: Any, allow_null_if_not_optional: bool) -> bool:
        if meta is None:
            return self.optional or allow_null_if_not_optional

        if isinstance(meta, str):
            return meta.lower() == "true" or meta.lower() == "false"
        if isinstance(meta, int):
            return meta == 0 or meta == 1

        return isinstance(meta, bool)

    def coerce_value_if_possible(self, value):
        if isinstance(value, bool):
            return value
        if isinstance(value, str):
            return value.lower() == "true"
        if isinstance(value, int):
            return value != 0

        raise Exception(f"Value {value} could not be coerced to boolean type")

    def invalid_value_hint(self, meta):
        if meta is None:
            return "value was null"
        if self.validate_value(meta, True):
            return None
        return f"Value was of type {type(meta)}, expected bool"


class File(DataType):
    def __init__(
        self, optional=False, extension=None, alternate_extensions: Set[str] = None
    ):
        """
        :param optional:
        :param extension: Used in CWL to try and guess the file extension where it's not available otherwise
        """
        super(File, self).__init__(optional=optional)
        self.extension = extension
        self.alternate_extensions = alternate_extensions

    def get_extensions(self):
        exts = []
        if self.extension:
            exts.append(self.extension)
        if self.alternate_extensions:
            exts.extend(self.alternate_extensions)

        return exts

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
        return {"class": "File", "path": value}

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

    def can_receive_from(self, other, source_has_default=False) -> bool:
        o = get_instantiated_type(other).received_type()
        if type(self) == File and isinstance(o, File):
            return True
        return super().can_receive_from(o)

    # def cwl_type(self, has_default=False):
    #     secs = self.secondary_files()
    #     if secs:
    #         tp = cwlgen.File(secondaryFiles=self.secondary_files())
    #         return [tp, "null"] if self.optional and not has_default else tp
    #     return super().cwl_type(has_default=has_default)


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
        return {"class": "Directory", "path": value}

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

    def is_array(self):
        return True

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
            type="array"
            # label=None,
            # input_binding=None
        )
        return [inp, "null"] if self.optional and not has_default else inp

    def map_cwl_type(self, parameter: cwlgen.Parameter) -> cwlgen.Parameter:
        parameter.type = cwlgen.CommandInputArraySchema(items=None, type="array")
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
        if other.is_array():
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
        if st.is_array():
            return st.fundamental_type()
        return st.received_type()

    def received_type(self):
        return Array(self._t.received_type(), optional=self.optional)


class Stdout(File):
    @staticmethod
    def name():
        return "Stdout"

    def __init__(self, subtype=None, optional=None):
        super().__init__(optional=False)

        subtype = get_instantiated_type(subtype) if subtype is not None else File()
        if optional is not None:
            subtype.optional = optional

        if subtype and not isinstance(subtype, File):
            raise Exception(
                "Janis does not currently support non-File stdout annotations"
            )

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
        st = self.subtype
        if self.optional is not None:
            st.optional = self.optional
        return st

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

    def __init__(self, subtype=None, stderrname=None, optional=None):
        super().__init__(optional=False)

        subtype = get_instantiated_type(subtype) if subtype is not None else File()
        if optional is not None:
            subtype.optional = optional

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
        st = self.subtype
        if self.optional is not None:
            st.optional = self.optional
        return st

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
    Stderr,
    Array,
]


def get_from_python_type(dt, optional: bool = None, overrider=None):
    if dt is None:
        return Boolean(optional=True)

    bc = overrider or get_instantiated_type
    dtt = dt if type(dt) == type else None
    typedt = type(dt)

    try:
        if dtt == str or typedt == str:
            return String(optional=optional)
    except Exception as e:
        print(e)
    if dtt == bool or typedt == bool:
        return Boolean(optional=optional)
    if dtt == int or typedt == int:
        return Int(optional=optional)
    if dtt == float or typedt == float:
        return Float(optional=optional)

    if is_qualified_generic(dt):

        if str(dt).startswith("typing.List"):
            nt = bc(dt.__args__[0], overrider=bc)
            return Array(nt, optional=optional)

        elif str(dt).startswith("typing.Union"):
            subtypes = dt.__args__
            # Filter out None or NoneType
            try:
                new_subtypes = [
                    t for t in subtypes if t is not None and type(None) != t
                ]
            except Exception as e:
                Logger.critical(
                    f"Couldn't determine the appropriate internal types from {str(dt)}, failed with error: {str(e)}"
                )
                raise
            optional = len(subtypes) != len(new_subtypes)

            if len(new_subtypes) == 0:
                raise TypeError(
                    "Unsure how to parse generic: '{str(dt)}', please raise an issue if you think this is in error"
                )

            if len(new_subtypes) == 1:
                return get_instantiated_type(
                    new_subtypes[0], optional=optional, overrider=bc
                )

            nts = [bc(n, overrider=bc) for n in new_subtypes]
            return UnionType(*nts, optional=optional)

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


NumericType = UnionType(Int, Double, Float)
AnyType = UnionType(String, Boolean, Int, Double, Float, File, Directory)
