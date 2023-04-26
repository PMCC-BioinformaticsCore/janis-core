

from dataclasses import dataclass
from janis_core import ToolInput

from janis_core import translation_utils as utils
from janis_core.translation_utils import DTypeType

@dataclass 
class Attributes:
    prefix:     bool
    optional:   bool
    default:    bool
    array:      bool
    prefixeach: bool

def get_attributes(tinput: ToolInput) -> Attributes:
    return Attributes(
        prefix=is_prefix(tinput),
        optional=is_optional(tinput),
        default=is_default(tinput),
        array=is_array(tinput),
        prefixeach=is_prefixeach(tinput)
    )

def is_prefix(tinput: ToolInput) -> bool:
    if tinput.prefix is not None:
        return True
    return False

def is_optional(tinput: ToolInput) -> bool:
    if tinput.input_type.optional == True: 
        return True
    return False

def is_default(tinput: ToolInput) -> bool:
    if tinput.default is not None: 
        return True
    return False

def is_array(tinput: ToolInput) -> bool:
    dtt = utils.get_dtt(tinput.input_type)
    if dtt in [
        DTypeType.SECONDARY_ARRAY,
        DTypeType.SECONDARY,
        DTypeType.FILE_PAIR_ARRAY,
        DTypeType.FILE_PAIR,
        DTypeType.FILE_ARRAY,
        DTypeType.GENERIC_ARRAY,
    ]:
        return True
    return False

def is_prefixeach(tinput: ToolInput) -> bool:
    if tinput.prefix_applies_to_all_elements == True:
        return True
    return False

