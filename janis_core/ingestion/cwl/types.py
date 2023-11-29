from typing import Any, Optional, Type
import regex as re
import sys
import inspect

from janis_core.ingestion.cwl.identifiers import get_id_entity
from janis_core.messages import ErrorCategory
from janis_core.messages import log_warning
from janis_core.utils.errors import UnsupportedError
from janis_core import settings

from janis_core import translation_utils as utils
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
from janis_core.redefinitions.types import (
    Fasta,
    Fastq,
    Sam,
    Bam,
    Bed,
    Vcf,
    Tsv,
    Cram,
)
from janis_core.ingestion.galaxy.datatypes.galaxy import (
    Anndata,
    BigWig,
    Gtf,
    Gff,
    BedGraph,
    BigBed
)

edam_format_map = {
    'edam:format_1929': Fasta,
    'edam:format_1930': Fastq,
    'edam:format_1931': Fastq,
    'edam:format_1932': Fastq,
    'edam:format_1933': Fastq,
    'edam:format_1974': Gff,
    'edam:format_2306': Gtf,
    'edam:format_2572': Bam,
    'edam:format_2573': Sam,
    'edam:format_3003': Bed,
    'edam:format_3004': BigBed,
    'edam:format_3006': BigWig,
    'edam:format_3016': Vcf,
    'edam:format_3462': Cram,
    'edam:format_3583': BedGraph,
    'edam:format_3590': Anndata,
    'edam:format_3709': Tsv,
}

file_priorities = {
    'Fasta': 1,
    'Fastq': 2,
    'Gff': 3,
    'Gtf': 4,
    'Bam': 5,
    'Sam': 6,
    'Bed': 7,
    'BigBed': 8,
    'BigWig': 9,
    'Vcf': 10,
    'Cram': 11,
    'BedGraph': 12,
    'Anndata': 13,
    'Tsv': 14,
}

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
    items = [x.replace('^', '') for x in items]
    return hash("|".join(sorted(set(items))))

def is_javascript_expr(text: str) -> bool:
    pattern = r'^\$((\([\s\S]*\))|(\{[\s\S]*\}))$'
    if re.search(pattern, text):
        return True
    return False


def ingest_cwl_type(
    cwl_type: Any, 
    cwl_utils: Any,
    cwl_entity: Any, 
    tool_uuid: str,
    secondaries: Optional[str | list[str]]=None,
    ) -> DataType:
    dtype_parser = CWLTypeParser(cwl_type, cwl_utils, cwl_entity, secondaries, tool_uuid)
    return dtype_parser.parse()


class CWLTypeParser:
    def __init__(
        self,
        cwl_type: Any,
        cwl_utils: Any,
        cwl_entity: Any, 
        secondaries: Optional[str | list[str]],
        tool_uuid: str
    ):
        self.cwl_type = cwl_type
        self.cwl_utils = cwl_utils
        self.cwl_entity = cwl_entity
        self.secondaries = secondaries
        self.tool_uuid = tool_uuid
        self.secondary_patterns: list[str] = self.preprocess_secondary_file_patterns(secondaries)

    @property 
    def entity_name(self) -> str:
        return get_id_entity(self.cwl_entity.id)
    
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
    def file_subtypes(self) -> list[Type[DataType]]:
        # TODO check if any format is a javascript expression.

        # get formats from cwl entity
        if not hasattr(self.cwl_entity, 'format') or self.cwl_entity.format is None:
            formats_list = []
        elif isinstance(self.cwl_entity.format, list):
            formats_list = self.cwl_entity.format
        elif isinstance(self.cwl_entity.format, str):
            formats_list = [self.cwl_entity.format]
        else:
            raise NotImplementedError

        # replace uri with edam:format_xxxx structure if applicable
        formats_list = [x.replace('http://edamontology.org/', 'edam:') for x in formats_list]
        # pull out edam formats
        edam_formats = [x for x in formats_list if x.startswith('edam:format_')]
        # remove unknown formats
        edam_formats = [x for x in edam_formats if x in edam_format_map]
        # return janis datatypes
        return [edam_format_map[x] for x in edam_formats]

    @property
    def secondary_files_plain(self) -> list[str]:
        return [x for x in self.secondary_patterns if not is_javascript_expr(x)]
    
    @property
    def secondary_files_exprs(self) -> list[str]:
        return [x for x in self.secondary_patterns if is_javascript_expr(x)]
    
    def parse(self) -> DataType:
        # parse basic cwl type
        inp_type = self.from_cwl_inner_type(self.cwl_type)
        return inp_type
    
    @property
    def file_subtype(self) -> DataType:
        filetypes = self.file_subtypes
        
        if len(filetypes) > 1:
            filetypes.sort(key=lambda x: file_priorities.get(x.__name__, 999))
            ftype = filetypes[0]
            msg = f'{self.entity_name}: supports multiple edam formats. selected {ftype.__name__} as fallback.'
            log_warning(self.tool_uuid, msg, ErrorCategory.DATATYPE)
            return ftype()
        
        elif len(filetypes) == 1:
            ftype = filetypes[0]
            return ftype()
        
        else:
            return File()
        
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
                    inner = self.file_subtype
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
                    msg = f"{self.entity_name}: unsupported datatype {cwl_type}. treated as generic File."
                    log_warning(self.tool_uuid, msg, ErrorCategory.DATATYPE)
                    inner = File()
                else:
                    raise UnsupportedError(f"can't detect type {cwl_type}")
            
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
                dtype = ingest_cwl_type(ctype, self.cwl_utils, self.cwl_entity, self.tool_uuid, self.secondaries)
                dtypes.append(dtype)
            
            # annotate janis dtypes as optional
            for dtype in dtypes:
                dtype.optional = is_optional
            
            # single type
            if len(dtypes) == 1:
                return dtypes[0]
            
            # multiple types
            elif len(dtypes) > 1:
                dtype_names = [x.name() for x in dtypes]
                edam_types = [x for x in dtypes if x.__class__.__name__ in edam_format_map]
                edam_types.sort(key=lambda x: file_priorities.get(x.__name__, 999))
                
                # edam file types
                if len(edam_types) > 0:
                    selected = edam_types[0]
                # other types
                else:
                    # create (dtype, dtypetype) tuples to label with type categories
                    dtypes_w_categories = [(utils.get_dtt(x), x) for x in dtypes]
                    # sort by order of importance (according to dtypetype category - uses enum index) 
                    dtypes_w_categories.sort(key=lambda x: x[0].value)
                    # get the most important type
                    selected = dtypes_w_categories[0][1]
                
                # print user message
                msg = f'{self.entity_name}: supports multiple datatypes {dtype_names}. selected {selected.__class__.__name__} as fallback.'
                log_warning(self.tool_uuid, msg, ErrorCategory.DATATYPE)
                return selected
            
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
                msg = f"{self.entity_name}: unsupported datatype {type(cwl_type).__name__}. Treated as a file."
                log_warning(self.tool_uuid, msg, ErrorCategory.DATATYPE)
                return File(optional=False)
            else:
                raise UnsupportedError(f"Can't parse type {type(cwl_type).__name__}")

    def get_data_type_from_secondaries(self) -> DataType:
        # TODO needs work here - add galaxy secondary types?
        sec_types = inspect.getmembers(sys.modules['janis_core.redefinitions.types'], inspect.isclass)
        sec_types = [x[1] for x in sec_types if issubclass(x[1], DataType)]
        sec_types = [x for x in sec_types if x.secondary_files()]
        sec_types_map = {}
        for dtype in sec_types:
            secondaries = dtype.secondary_files()
            assert(isinstance(secondaries, list))
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
            msg = f'{self.entity_name}: Could not parse datatype from javascript expression. Treated as generic File with secondaries.'
            log_warning(self.tool_uuid, msg, ErrorCategory.DATATYPE)
            return GenericFileWithSecondaries(secondaries=[])
        
        else:
            raise NotImplementedError

    
            

    

