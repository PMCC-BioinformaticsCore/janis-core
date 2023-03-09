

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
    from ruamel.yaml.scalarstring import DoubleQuotedScalarString

    if isinstance(cwlvalue, DoubleQuotedScalarString):
        return str(cwlvalue)
    elif isinstance(cwlvalue, CommentedSeq):
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
        self.patterns = self.preprocess_secondary_file_patterns(secondary_files)
    
    def preprocess_secondary_file_patterns(self, secondary_files: Optional[list[str]]) -> Optional[list[str]]:
        patterns: Optional[list[str]] = None
        # preprocessing to get patterns for secondary files
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
        if self.patterns:
            if len(self.patterns) > 0 and not self.patterns[0].startswith('$('):
                return self.patterns
        return None
    
    @property
    def secondary_files_expr(self) -> Optional[str]:
        if self.patterns:
            if len(self.patterns) == 1 and self.patterns[0].startswith('$('):
                return self.patterns[0]
        return None
    
    def parse(self) -> Tuple[DataType, list[str]]:
        inp_type = self.from_cwl_inner_type(self.cwl_type)
        
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
                inp_type = GenericFileWithSecondaries(secondaries=[])

        return (inp_type, self.error_msgs)
    
    def from_cwl_inner_type(self, cwl_type: Any) -> DataType:
        if isinstance(cwl_type, str):
            optional = "?" in cwl_type
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
            return inner(optional=optional)

        elif isinstance(cwl_type, list):
            optional = None
            types = []
            for c in cwl_type:
                if c == "null":
                    optional = True
                else:
                    dtype, error_messages = ingest_cwl_type(c, self.cwl_utils, [])
                    self.error_msgs += error_messages
                    types.append(dtype)

            if len(types) == 1:
                if optional is not None:
                    types[0].optional = optional
                return types[0]
            else:
                from janis_core.types.common_data_types import UnionType

                if optional is not None:
                    for inner in types:
                        inner.optional = optional

                return UnionType(*types)

        elif isinstance(cwl_type, self.cwl_utils.CommandInputArraySchema):
            return Array(self.from_cwl_inner_type(cwl_type.items))
        elif isinstance(cwl_type, self.cwl_utils.InputArraySchema):
            return Array(self.from_cwl_inner_type(cwl_type.items))
        elif isinstance(cwl_type, self.cwl_utils.CommandOutputArraySchema):
            return Array(self.from_cwl_inner_type(cwl_type.items))
        elif isinstance(cwl_type, self.cwl_utils.OutputArraySchema):
            return Array(self.from_cwl_inner_type(cwl_type.items))
        elif isinstance(cwl_type, self.cwl_utils.InputEnumSchema):
            return String()

        else:
            if settings.datatypes.ALLOW_UNPARSEABLE_DATATYPES:
                msg = f"Unsupported datatype: {type(cwl_type).__name__}. Treated as a file."
                self.error_msgs.append(msg)
                return File(optional=False)
            else:
                raise UnsupportedError(f"Can't parse type {type(cwl_type).__name__}")

    def get_data_type_from_secondaries(self, optional: bool) -> DataType:
        global file_datatype_cache

        if not file_datatype_cache:
            FastaGzType = None
            try:
                from janis_bioinformatics.data_types import FastaGz

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
            return file_datatype_cache[sec_hash](optional=optional)

        return GenericFileWithSecondaries(secondaries=secondary_files)

    
            

    

