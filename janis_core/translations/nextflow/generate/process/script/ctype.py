

from enum import Enum, auto

from janis_core import ToolInput
from janis_core.types import Boolean, Filename, Array


# public enum
class CType(Enum):

    FLAG_BASIC           = auto()
    FLAG_DEFAULT         = auto()
    FLAG_OPTIONAL        = auto()

    FLAG_ARR_BASIC           = auto()
    FLAG_ARR_DEFAULT         = auto()
    FLAG_ARR_OPTIONAL        = auto()

    POS_BASIC           = auto()
    POS_DEFAULT         = auto()
    POS_OPTIONAL        = auto()

    POS_BASIC_ARR       = auto()
    POS_DEFAULT_ARR     = auto()
    POS_OPTIONAL_ARR    = auto()

    POS_BASIC_ARR_PREFIXEACH    = auto()
    POS_DEFAULT_ARR_PREFIXEACH  = auto()
    POS_OPTIONAL_ARR_PREFIXEACH = auto()

    OPT_BASIC           = auto()
    OPT_DEFAULT         = auto()
    OPT_OPTIONAL        = auto()

    OPT_BASIC_ARR               = auto()
    OPT_DEFAULT_ARR             = auto()
    OPT_OPTIONAL_ARR            = auto()
    
    OPT_BASIC_ARR_PREFIXEACH    = auto()
    OPT_DEFAULT_ARR_PREFIXEACH  = auto()
    OPT_OPTIONAL_ARR_PREFIXEACH = auto()


# public function
def get_ctype(tinput: ToolInput) -> CType:

    # flag arrays
    if is_array_type(tinput) and is_flag(tinput) and is_basic(tinput):
        return CType.FLAG_ARR_BASIC
    elif is_array_type(tinput) and is_flag(tinput) and has_default(tinput):
        return CType.FLAG_ARR_DEFAULT
    elif is_array_type(tinput) and is_flag(tinput) and is_optional(tinput):
        return CType.FLAG_ARR_OPTIONAL
    
    # flags
    elif is_flag(tinput) and is_basic(tinput):
        return CType.FLAG_ARR_BASIC
    elif is_flag(tinput) and has_default(tinput):
        return CType.FLAG_ARR_DEFAULT
    elif is_flag(tinput) and is_optional(tinput):
        return CType.FLAG_ARR_OPTIONAL
    
    # positional arrays
    elif is_array_type(tinput) and is_positional(tinput) and is_basic(tinput):
        return CType.POS_BASIC_ARR
    elif is_array_type(tinput) and is_positional(tinput) and has_default(tinput):
        return CType.POS_DEFAULT_ARR
    elif is_array_type(tinput) and is_positional(tinput) and is_optional(tinput):
        return CType.POS_OPTIONAL_ARR
    
    # option arrays
    elif is_array_type(tinput) and is_option(tinput) and is_basic(tinput) and tinput.prefix_applies_to_all_elements:
        return CType.OPT_BASIC_ARR_PREFIXEACH
    elif is_array_type(tinput) and is_option(tinput) and has_default(tinput) and tinput.prefix_applies_to_all_elements:
        return CType.OPT_DEFAULT_ARR_PREFIXEACH
    elif is_array_type(tinput) and is_option(tinput) and is_optional(tinput) and tinput.prefix_applies_to_all_elements:
        return CType.OPT_OPTIONAL_ARR_PREFIXEACH
    elif is_array_type(tinput) and is_option(tinput) and is_basic(tinput):
        return CType.OPT_BASIC_ARR
    elif is_array_type(tinput) and is_option(tinput) and has_default(tinput):
        return CType.OPT_DEFAULT_ARR
    elif is_array_type(tinput) and is_option(tinput) and is_optional(tinput):
        return CType.OPT_OPTIONAL_ARR
    
    # positionals
    elif is_positional(tinput) and is_basic(tinput):
        return CType.POS_BASIC
    elif is_positional(tinput) and has_default(tinput):
        return CType.POS_DEFAULT
    elif is_positional(tinput) and is_optional(tinput):
        return CType.POS_OPTIONAL
    
    # options
    elif is_option(tinput) and is_basic(tinput):
        return CType.OPT_BASIC
    elif is_option(tinput) and has_default(tinput):
        return CType.OPT_DEFAULT
    elif is_option(tinput) and is_optional(tinput):
        return CType.OPT_OPTIONAL

    else:
        # some future CType
        raise NotImplementedError


### HELPER FUNCS ###

def is_positional(tinput: ToolInput) -> bool:
    if tinput.prefix == None:
        return True
    return False

def is_flag(tinput: ToolInput) -> bool:
    if isinstance(tinput.input_type, Boolean):
        return True
    return False

def is_option(tinput: ToolInput) -> bool:
    if not is_positional(tinput) and not is_flag(tinput):
        return True
    return False

def is_basic(tinput: ToolInput) -> bool:
    if not is_optional(tinput):
        if tinput.default is None:
            return True
    return False

def is_optional(tinput: ToolInput) -> bool:
    if isinstance(tinput.input_type, Filename):
        return False
    elif tinput.input_type.optional == True:  # type: ignore
        return True
    return False

def has_default(tinput: ToolInput) -> bool:
    if tinput.default is not None:
        return True
    return False

def is_flag_true(tinput: ToolInput) -> bool:
    if is_flag(tinput) and str(tinput.default) != 'False':
        return True
    return False

def is_flag_false(tinput: ToolInput) -> bool:
    if is_flag(tinput) and not is_flag_true(tinput):
        return True
    return False

def is_array_type(tinput: ToolInput) -> bool:
    if isinstance(tinput.input_type, Array) and tinput.input_type.name() == 'Array':
        return True
    return False

    