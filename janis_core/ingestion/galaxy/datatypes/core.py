
# oi scroll to bottom 

from .JanisDatatype import JanisDatatype

file_t = JanisDatatype(
    format='file',
    source='janis',
    classname='File',
    extensions=None,
    import_path='janis_core.types.common_data_types'
)
string_t = JanisDatatype(
    format='string',
    source='janis',
    classname='String',
    extensions=None,
    import_path='janis_core.types.common_data_types'
)
directory_t = JanisDatatype(
    format='directory',
    source='janis',
    classname='Directory',
    extensions=None,
    import_path='janis_core.types.common_data_types'
)
int_t = JanisDatatype(
    format='integer',
    source='janis',
    classname='Int',
    extensions=None,
    import_path='janis_core.types.common_data_types'
)
float_t = JanisDatatype(
    format='float',
    source='janis',
    classname='Float',
    extensions=None,
    import_path='janis_core.types.common_data_types'
)
bool_t = JanisDatatype(
    format='boolean',
    source='janis',
    classname='Boolean',
    extensions=None,
    import_path='janis_core.types.common_data_types'
)

DEFAULT_DATATYPE = file_t

CORE_DATATYPES = set([
    file_t.classname,
    string_t.classname,
    directory_t.classname,
    int_t.classname,
    float_t.classname,
    bool_t.classname
])

