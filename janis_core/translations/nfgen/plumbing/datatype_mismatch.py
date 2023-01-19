

from janis_core.types import DataType

from .. import nfgen_utils


# public
# check if mismatch is occuring between source and dest 
def is_datatype_mismatch(srctype: DataType, desttype: DataType) -> bool:
    if _array_type(srctype) and _single_type(desttype):
        return True
    elif _single_type(srctype) and _array_type(desttype):
        return True
    elif _secondary_array_type(srctype) and _secondary_type(desttype):
        return True
    elif _secondary_type(srctype) and _secondary_array_type(desttype):
        return True
    return False


# public
# handle a mismatch we have encountered
# each function returns an expression we can tack onto channel 
# (during process / subworkflow call) to fix
def handle_datatype_mismatch(srctype: DataType, desttype: DataType) -> str:
    "delegate the mismatch to specific handling function"
    if _array_type(srctype) and _single_type(desttype):
        return '.flatten().first()'
    
    elif _single_type(srctype) and _array_type(desttype):
        return '.toList()'
    
    elif _secondary_array_type(srctype) and _secondary_type(desttype):
        basetype = nfgen_utils.get_base_type(srctype)
        exts = nfgen_utils.get_extensions(basetype)
        size = len(exts)
        return f'.flatten().collate( {size} ).first()'
    
    elif _secondary_type(srctype) and _secondary_array_type(desttype):
        return '.flatten().collect()'
    
    else:
        raise RuntimeError('DEV: there should be a type mismatch here, but apparently not?')


# private 
# helper functions to identify types involved in mismatches
def _single_type(dtype: DataType) -> bool:
    if not dtype.is_array():
        if not _secondary_type(dtype):
            return True
    return False

def _array_type(dtype: DataType) -> bool:
    if not nfgen_utils.is_array_secondary_type(dtype):
        if dtype.is_array():
            return True
    return False

def _secondary_type(dtype: DataType) -> bool:
    return nfgen_utils.is_secondary_type(dtype)

def _secondary_array_type(dtype: DataType) -> bool:
    return nfgen_utils.is_array_secondary_type(dtype)