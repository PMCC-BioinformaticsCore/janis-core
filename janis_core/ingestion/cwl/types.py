

from typing import Any, Optional, Tuple
import regex as re
from copy import deepcopy
import sys
import inspect

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

def _calcluate_hash_of_set(items: Any):
    return hash("|".join(sorted(set(items))))


def is_javascript_expr(text: str) -> bool:
    pattern = r'^\$((\([\s\S]*\))|(\{[\s\S]*\}))$'
    if re.search(pattern, text):
        return True
    return False


def ingest_cwl_type(
    cwl_type: Any, 
    cwl_utils: Any,
    secondaries: Optional[str | list[str]]=None,
    ) -> Tuple[DataType, list[str]]:
    dtype_parser = CWLTypeParser(cwl_type, cwl_utils, secondaries)
    return dtype_parser.parse()


class CWLTypeParser:
    def __init__(
        self,
        cwl_type: Any,
        cwl_utils: Any,
        secondaries: Optional[str | list[str]],
    ):
        self.cwl_type = cwl_type
        self.cwl_utils = cwl_utils
        self.secondaries = secondaries
        self.error_msgs: list[str] = []
        self.secondary_patterns: list[str] = self.preprocess_secondary_file_patterns(secondaries)
    
    def preprocess_secondary_file_patterns(self, secondaries: Optional[str | list[str]]) -> list[str]:
        out: list[str] = []

        if not secondaries:
            return out
        
        elif isinstance(secondaries, str):
            out.append(secondaries.strip())
            return out
        
        for sfile in secondaries:
            if hasattr(sfile, 'pattern'):
                out.append(sfile.pattern.strip())
            else:
                out.append(sfile.strip())

        return out

    @property
    def secondary_files_plain(self) -> list[str]:
        return [x for x in self.secondary_patterns if not is_javascript_expr(x)]
    
    @property
    def secondary_files_exprs(self) -> list[str]:
        return [x for x in self.secondary_patterns if is_javascript_expr(x)]
    
    def parse(self) -> Tuple[DataType, list[str]]:
        # parse basic cwl type
        inp_type = self.from_cwl_inner_type(self.cwl_type)
        return (inp_type, self.error_msgs)
    
    def from_cwl_inner_type(self, cwl_type: Any) -> DataType:
        if isinstance(cwl_type, str):
            # optionality
            is_optional = "?" in cwl_type
            cwl_type = cwl_type.replace("?", "")
            
            # array depth - is this even right?
            array_count = 0
            while cwl_type.endswith("[]"):  # TODO is this the right syntax? [][] vs [[]]
                array_count += 1
                cwl_type = cwl_type[:-2]

            # inner type
            if cwl_type == "File":
                if self.secondary_patterns:
                    inner = self.get_data_type_from_secondaries()
                else:
                    inner = File()
            elif cwl_type == "Directory":
                inner = Directory()
            elif cwl_type == "string":
                inner = String()
            elif cwl_type == "int":
                inner = Int()
            elif cwl_type == "float":
                inner = Float()
            elif cwl_type == "double":
                inner = Float()
            elif cwl_type == "boolean":
                inner = Boolean()
            elif cwl_type == "stdout":
                inner = Stdout()
            elif cwl_type == "stderr":
                inner = Stderr()
            elif cwl_type == "Any":
                inner = String()
            elif cwl_type == "long":
                inner = Int()
            else:
                if settings.datatypes.ALLOW_UNPARSEABLE_DATATYPES:
                    msg = f"Unsupported datatype: {cwl_type}. Treated as a file."
                    self.error_msgs.append(msg)
                    inner = File()
                else:
                    raise UnsupportedError(f"Can't detect type {cwl_type}")
            
            # array depth
            dtype = inner
            for i in range(array_count):
                dtype = Array(dtype)
            dtype.optional = is_optional
            return dtype

        elif isinstance(cwl_type, list):
            # optionality
            is_optional = True if 'null' in cwl_type else False
            
            # individual cwl types
            cwl_types = [x for x in cwl_type if x != 'null']
            
            # casting individual cwl types to janis
            dtypes: list[DataType] = []
            for ctype in cwl_types:
                dtype, error_messages = ingest_cwl_type(ctype, self.cwl_utils, self.secondaries)
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

    def get_data_type_from_secondaries(self) -> DataType:
        # TODO needs work here - add galaxy secondary types
        sec_types = inspect.getmembers(sys.modules['janis_core.redefinitions.types'], inspect.isclass)
        sec_types = [x[1] for x in sec_types if issubclass(x[1], DataType)]
        sec_types = [x for x in sec_types if x.secondary_files()]
        sec_types_map = {}
        for dtype in sec_types:
            secondaries = dtype.secondary_files()
            assert(isinstance(secondaries, list))
            secondaries = [x.replace('^', '') for x in secondaries]
            sec_hash = _calcluate_hash_of_set(secondaries)
            sec_types_map[sec_hash] = dtype

        if self.secondary_files_plain and self.secondary_files_exprs:
            raise NotImplementedError
        
        elif self.secondary_files_plain:
            sec_hash = _calcluate_hash_of_set(self.secondary_files_plain)
            if sec_hash in sec_types_map:
                dtype = sec_types_map[sec_hash]
                return dtype()
            else:
                return GenericFileWithSecondaries(secondaries=self.secondary_files_plain)
        
        elif self.secondary_files_exprs:
            # res, success = parse_basic_expression(self.secondary_files_exprs)
            # if success:
            #     raise NotImplementedError
            # else:
            msg = f'could not parse secondaries format from javascript expressions: {self.secondary_files_exprs[0]}'
            self.error_msgs.append(msg)
            return GenericFileWithSecondaries(secondaries=[])
        
        else:
            raise NotImplementedError

    
            

    

