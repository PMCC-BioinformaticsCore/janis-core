

from enum import Enum, auto

from janis_core import ToolInput
from janis_core.types import Boolean, Filename

"""
This file labels a ToolInput as a particular IType using get_itype(). 
This allows us to know how the tool input should be formatted in the prescript
and script section of the process. 
"""


"""
POS_BASIC
A basic positional. Has no prefix, and is mandatory.
Will have either a process input or param input, or will be fed a value via InputSelector.
"""

class IType(Enum):
    FLAG_TRUE           = auto()
    FLAG_FALSE          = auto()

    POS_BASIC           = auto()
    POS_DEFAULT         = auto()
    POS_OPTIONAL        = auto()

    POS_BASIC_ARR       = auto()
    POS_DEFAULT_ARR     = auto()
    POS_OPTIONAL_ARR    = auto()

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
def get_itype(tinput: ToolInput) -> IType:
    # flags
    if is_flag_true(tinput):
        return IType.FLAG_TRUE
    elif is_flag_false(tinput):
        return IType.FLAG_FALSE

    # positional arrays
    elif is_array(tinput) and is_positional(tinput) and is_basic(tinput):
        return IType.POS_BASIC_ARR
    elif is_array(tinput) and is_positional(tinput) and has_default(tinput):
        return IType.POS_DEFAULT_ARR
    elif is_array(tinput) and is_positional(tinput) and is_optional(tinput):
        return IType.POS_OPTIONAL_ARR
    
    # option arrays
    elif is_array(tinput) and is_option(tinput) and is_basic(tinput) and tinput.prefix_applies_to_all_elements:
        return IType.OPT_BASIC_ARR_PREFIXEACH
    elif is_array(tinput) and is_option(tinput) and has_default(tinput) and tinput.prefix_applies_to_all_elements:
        return IType.OPT_DEFAULT_ARR_PREFIXEACH
    elif is_array(tinput) and is_option(tinput) and is_basic(tinput) and tinput.prefix_applies_to_all_elements:
        return IType.OPT_OPTIONAL_ARR_PREFIXEACH
    elif is_array(tinput) and is_option(tinput) and is_basic(tinput):
        return IType.OPT_BASIC_ARR
    elif is_array(tinput) and is_option(tinput) and has_default(tinput):
        return IType.OPT_DEFAULT_ARR
    elif is_array(tinput) and is_option(tinput) and is_optional(tinput):
        return IType.OPT_OPTIONAL_ARR
    
    # positionals
    elif is_positional(tinput) and is_basic(tinput):
        return IType.POS_BASIC
    elif is_positional(tinput) and has_default(tinput):
        return IType.POS_DEFAULT
    elif is_positional(tinput) and is_optional(tinput):
        return IType.POS_OPTIONAL
    
    # options
    elif is_option(tinput) and is_basic(tinput):
        return IType.OPT_BASIC
    elif is_option(tinput) and has_default(tinput):
        return IType.OPT_DEFAULT
    elif is_option(tinput) and is_optional(tinput):
        return IType.OPT_OPTIONAL

    else:
        # some future IType
        raise NotImplementedError



# bool ToolInput type identities
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


# bool features
def is_basic(tinput: ToolInput) -> bool:
    if not is_optional(tinput):
        if tinput.default == None:
            return True
    return False

def is_optional(tinput: ToolInput) -> bool:
    if isinstance(tinput.input_type, Filename):
        return False
    elif tinput.input_type.optional == True:
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

def is_array(tinput: ToolInput) -> bool:
    if tinput.input_type.is_array():
        return True
    return False
