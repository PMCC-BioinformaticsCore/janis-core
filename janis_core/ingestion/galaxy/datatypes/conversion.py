

from .JanisDatatype import JanisDatatype
from .core import bool_t
from .core import file_t
from .core import string_t
from .core import directory_t
from .core import float_t
from .core import int_t
from .core import CORE_DATATYPES
from .core import DEFAULT_DATATYPE
from .register import register

core_type_priorities = {
    bool_t.classname: 1,  # highest priority
    directory_t.classname: 2,
    float_t.classname: 3,
    string_t.classname: 4,
    int_t.classname: 5,
    file_t.classname: 6, # lowest priority
}

def galaxy_to_janis(galaxy_types: list[str]) -> list[JanisDatatype]:
    out: list[JanisDatatype] = []
    for gtype in galaxy_types:
        jtype = register.get(gtype)
        if jtype is not None:
            out.append(jtype)
    return out

def janis_to_core(query_types: list[JanisDatatype]) -> list[JanisDatatype]:
    core: dict[str, JanisDatatype] = {} # dict to keep unique
    for qtype in query_types:
        if qtype.classname in CORE_DATATYPES:
            core[qtype.classname] = qtype
    return list(core.values())

def select_primary_core_type(query_types: list[JanisDatatype]) -> JanisDatatype:
    if query_types:
        selected_type = query_types[0]
        for jtype in query_types[1:]:
            if core_type_priorities[jtype.classname] > core_type_priorities[selected_type.classname]:
                selected_type = jtype
        return selected_type
    else:
        return DEFAULT_DATATYPE
