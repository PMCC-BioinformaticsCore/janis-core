

from janis_core.types import DataType

from .. import nfgen_utils

from .common import array_type
from .common import single_type
from .common import secondary_type
from .common import secondary_array_type


# public
# check if mismatch is occuring between source and dest 
def is_datatype_mismatch(srctype: DataType, desttype: DataType) -> bool:
    if array_type(srctype) and single_type(desttype):
        return True
    elif single_type(srctype) and array_type(desttype):
        return True
    elif secondary_array_type(srctype) and secondary_type(desttype):
        return True
    elif secondary_type(srctype) and secondary_array_type(desttype):
        return True
    return False


# public
# handle a mismatch we have encountered
# each function returns an expression we can tack onto channel 
# (during process / subworkflow call) to fix
def handle_datatype_mismatch(srctype: DataType, desttype: DataType) -> str:
    "delegate the mismatch to specific handling function"
    if array_type(srctype) and single_type(desttype):
        return '.flatten().first()'
    
    elif single_type(srctype) and array_type(desttype):
        return '.toList()'
    
    elif secondary_array_type(srctype) and secondary_type(desttype):
        basetype = nfgen_utils.get_base_type(srctype)
        exts = nfgen_utils.get_extensions(basetype)
        size = len(exts)
        return f'.flatten().collate( {size} ).first()'
    
    elif secondary_type(srctype) and secondary_array_type(desttype):
        return '.flatten().toList()'
    
    else:
        raise RuntimeError('DEV: there should be a type mismatch here, but apparently not?')



