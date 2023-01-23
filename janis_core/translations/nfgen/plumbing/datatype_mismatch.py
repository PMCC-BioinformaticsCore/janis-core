

from janis_core.types import DataType

from .common import array_type
from .common import single_type
from .common import secondary_type
from .common import secondary_array_type
from .common import get_collate_size


# public
# check if mismatch is occuring between source and dest 
def requires_data_operation(srctype: DataType, desttype: DataType, src_scatter: bool, dest_scatter: bool) -> bool:
    if secondary_array_type(desttype):
        return True
    elif secondary_array_type(srctype) and secondary_type(desttype):
        return True
    if array_type(srctype) and single_type(desttype):
        return True
    elif single_type(srctype) and array_type(desttype):
        return True
    return False


# public
# handle a mismatch we have encountered
# each function returns an expression we can tack onto channel 
# (during process / subworkflow call) to fix
def handle_data_operation(srctype: DataType, desttype: DataType, src_scatter: bool, dest_scatter: bool) -> str:
    # handle secondary array process input format
    if secondary_array_type(desttype):
        operation = '.flatten().toList()'
    
    elif array_type(srctype) and single_type(desttype):
        operation = '.flatten().first()'
    
    elif single_type(srctype) and array_type(desttype):
        operation = '.toList()'
    
    elif secondary_array_type(srctype) and secondary_type(desttype):
        size = get_collate_size(srctype)
        operation = f'.flatten().collate( {size} ).first()'
    
    else:
        raise RuntimeError('DEV: there should be a type mismatch here, but apparently not?')

    return operation
