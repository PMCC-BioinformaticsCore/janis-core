

from typing import Any, Optional, Tuple

from janis_core import JanisShed
from janis_core.types import (
    DataType, 
    GenericFileWithSecondaries, 
    File,
    Array,
    Directory,
    String,
    Int,
    Float,
    Boolean,
    Stdout,
    Stderr,
)

from janis_core.utils.errors import UnsupportedError
from janis_core import settings

from .expressions import parse_basic_expression


file_datatype_cache: dict[int, Any] = {}


def cast_cwl_type_to_python(cwlvalue: Any) -> Any:
    """
    cwl utils parser stores these as ruamel types, not python types.
    sometimes need to be converted to python for janis datatype inference.
    for example:
    often we see <class 'ruamel.yaml.comments.CommentedSeq'> for lists of strings
    each item being a <class 'ruamel.yaml.scalarstring.DoubleQuotedScalarString'>
    """
    from ruamel.yaml.comments import CommentedSeq
    from ruamel.yaml.comments import CommentedMap
    from ruamel.yaml.scalarstring import DoubleQuotedScalarString
    from ruamel.yaml.scalarstring import SingleQuotedScalarString
    from ruamel.yaml.scalarstring import ScalarString
    from ruamel.yaml.scalarint import ScalarInt
    from ruamel.yaml.scalarbool import ScalarBoolean
    from ruamel.yaml.scalarfloat import ScalarFloat

    if isinstance(cwlvalue, DoubleQuotedScalarString):
        return str(cwlvalue)
    elif isinstance(cwlvalue, SingleQuotedScalarString):
        return str(cwlvalue)
    elif isinstance(cwlvalue, ScalarString):
        return str(cwlvalue)
    elif isinstance(cwlvalue, ScalarInt):
        return int(cwlvalue)
    elif isinstance(cwlvalue, ScalarFloat):
        return float(cwlvalue)
    elif isinstance(cwlvalue, ScalarBoolean):
        return bool(cwlvalue)
    elif isinstance(cwlvalue, CommentedMap):
        out: dict[str, Any] = {}
        for key, val in cwlvalue.items():
            key = cast_cwl_type_to_python(key)
            val = cast_cwl_type_to_python(val)
            out[key] = val
        return out
    elif isinstance(cwlvalue, CommentedSeq):
        return [cast_cwl_type_to_python(x) for x in cwlvalue]
    elif isinstance(cwlvalue, list):
        return [cast_cwl_type_to_python(x) for x in cwlvalue]
    
    return cwlvalue

def _calcluate_hash_of_set(the_set: Any):
    return hash("|".join(sorted(set(the_set))))



def ingest_cwl_type(
    cwl_type: Any, 
    cwl_utils: Any,
    secondary_files: Optional[list[str]]=None,
    ) -> Tuple[DataType, list[str]]:
    dtype_parser = CWLTypeParser(cwl_type, cwl_utils, secondary_files)
    return dtype_parser.parse()


class CWLTypeParser:
    def __init__(
        self,
        cwl_type: Any,
        cwl_utils: Any,
        secondary_files: Optional[list[str]],
    ):
        self.cwl_type = cwl_type
        self.cwl_utils = cwl_utils
        self.error_msgs: list[str] = []
        self.secondary_patterns = self.preprocess_secondary_file_patterns(secondary_files)
    
    def preprocess_secondary_file_patterns(self, secondary_files: Optional[list[str]]) -> Optional[list[str]]:
        # preprocessing to get patterns for secondary files
        if isinstance(secondary_files, str):
            secondary_files = [secondary_files]

        patterns: Optional[list[str]] = None
        if secondary_files:
            patterns = []
            for sfile in secondary_files:
                if hasattr(sfile, 'pattern'):
                    patterns.append(sfile.pattern)
                else:
                    patterns.append(sfile)
        return patterns

    @property
    def secondary_files(self) -> Optional[list[str]]:
        if self.secondary_patterns:
            if len(self.secondary_patterns) > 0 and not self.secondary_patterns[0].startswith('$('):
                return self.secondary_patterns
        return None
    
    @property
    def secondary_files_expr(self) -> Optional[str]:
        if self.secondary_patterns:
            if len(self.secondary_patterns) == 1 and self.secondary_patterns[0].startswith('$('):
                return self.secondary_patterns[0]
        return None
    
    def parse(self) -> Tuple[DataType, list[str]]:
        inp_type = self.from_cwl_inner_type(self.cwl_type)
        
        # secondary files
        if self.secondary_files is not None:
            array_optional_layers: list[bool] = []
            while isinstance(inp_type, Array):
                array_optional_layers.append(inp_type.optional)
                inp_type = inp_type.subtype()

            inp_type = self.get_data_type_from_secondaries(inp_type.optional)
            for is_optional in array_optional_layers[::-1]:
                inp_type = Array(inp_type, optional=is_optional)
        
        elif self.secondary_files_expr is not None:
            res, success = parse_basic_expression(self.secondary_files_expr)
            if success:
                raise NotImplementedError
            else:
                msg = f'could not parse secondaries format from javascript expression: {res}'
                self.error_msgs.append(msg)
                inp_type = GenericFileWithSecondaries(secondaries=[], optional=inp_type.optional)

        return (inp_type, self.error_msgs)
    
    def from_cwl_inner_type(self, cwl_type: Any) -> DataType:
        if isinstance(cwl_type, str):
            is_optional = "?" in cwl_type
            cwl_type = cwl_type.replace("?", "")
            array_count = 0
            while cwl_type.endswith("[]"):
                array_count += 1
                cwl_type = cwl_type[:-2]

            if cwl_type == "File":
                inner = File
            elif cwl_type == "Directory":
                inner = Directory
            elif cwl_type == "string":
                inner = String
            elif cwl_type == "int":
                inner = Int
            elif cwl_type == "float":
                inner = Float
            elif cwl_type == "double":
                inner = Float
            elif cwl_type == "boolean":
                inner = Boolean
            elif cwl_type == "stdout":
                inner = Stdout
            elif cwl_type == "stderr":
                inner = Stderr
            elif cwl_type == "Any":
                inner = String
            elif cwl_type == "long":
                inner = Int
            else:
                if settings.datatypes.ALLOW_UNPARSEABLE_DATATYPES:
                    msg = f"Unsupported datatype: {cwl_type}. Treated as a file."
                    self.error_msgs.append(msg)
                    inner = File
                else:
                    raise UnsupportedError(f"Can't detect type {cwl_type}")
            return inner(optional=is_optional)

        elif isinstance(cwl_type, list):
            # optionality
            is_optional = True if 'null' in cwl_type else False
            
            # individual cwl types
            cwl_types = [x for x in cwl_type if x != 'null']
            
            # casting individual cwl types to janis
            dtypes: list[DataType] = []
            for ctype in cwl_types:
                dtype, error_messages = ingest_cwl_type(ctype, self.cwl_utils, []) # [] may be an mistake?
                self.error_msgs += error_messages
                dtypes.append(dtype)
            
            # annotate janis dtypes as optional
            for dtype in dtypes:
                dtype.optional = is_optional
            
            # single type
            if len(dtypes) == 1:
                return dtypes[0]
            
            # multiple types type
            elif len(dtypes) > 1:
                dtype_names = [x.name() for x in dtypes]
                msg = f'entity supports multiple datatypes: {dtype_names}. selected {dtype_names[0]} as fallback. this may affect pipeline execution'
                self.error_msgs.append(msg)
                return dtypes[0]
            
            else:
                raise RuntimeError

        elif isinstance(cwl_type, self.cwl_utils.CommandInputArraySchema):
            return Array(self.from_cwl_inner_type(cwl_type.items))
        elif isinstance(cwl_type, self.cwl_utils.InputArraySchema):
            return Array(self.from_cwl_inner_type(cwl_type.items))
        elif isinstance(cwl_type, self.cwl_utils.CommandOutputArraySchema):
            return Array(self.from_cwl_inner_type(cwl_type.items))
        elif isinstance(cwl_type, self.cwl_utils.OutputArraySchema):
            return Array(self.from_cwl_inner_type(cwl_type.items))
        elif isinstance(cwl_type, self.cwl_utils.InputEnumSchema):
            # NOTE does this require optionality checking?
            return String()

        else:
            if settings.datatypes.ALLOW_UNPARSEABLE_DATATYPES:
                msg = f"Unsupported datatype: {type(cwl_type).__name__}. Treated as a file."
                self.error_msgs.append(msg)
                return File(optional=False)
            else:
                raise UnsupportedError(f"Can't parse type {type(cwl_type).__name__}")

    def get_data_type_from_secondaries(self, is_optional: bool) -> DataType:
        global file_datatype_cache

        # TODO needs work here - JanisShed.get_all_datatypes() doesn't return anything

        if not file_datatype_cache:
            FastaGzType = None
            try:
                from janis_core.redefinitions.types import FastaGz
                # from janis_bioinformatics.data_types import FastaGz
                FastaGzType = FastaGz
            except ImportError:
                pass
            dts = JanisShed.get_all_datatypes()
            file_datatype_cache = {
                _calcluate_hash_of_set(dt.secondary_files()): dt
                for dt in dts
                if issubclass(dt, File)
                and dt().secondary_files()
                and (FastaGzType is not None and not issubclass(dt, FastaGzType))
            }

        secondary_files: list[str] = self.secondary_files # type: ignore
        sec_hash = _calcluate_hash_of_set(secondary_files)
        if sec_hash in file_datatype_cache:
            return file_datatype_cache[sec_hash](optional=is_optional)

        return GenericFileWithSecondaries(secondaries=secondary_files, optional=is_optional)

    
            

    

