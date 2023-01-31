

from janis_core.types import DataType

from .. import nfgen_utils

from .common import array_type
from .common import single_type
from .common import secondary_type
from .common import secondary_array_type
from .common import get_collate_size


# public
def is_datatype_mismatch(srctype: DataType, desttype: DataType) -> bool:
    if secondary_secondary_mismatch(srctype, desttype):
        return True
    elif secondary_single_mismatch(srctype, desttype):
        return True
    elif secondary_array_type(desttype):
        return True
    elif secondary_array_type(srctype) and secondary_type(desttype):
        return True
    elif array_type(srctype) and single_type(desttype):
        return True
    elif single_type(srctype) and array_type(desttype):
        return True
    return False

def generate_datatype_mismatch_plumbing(srctype: DataType, desttype: DataType) -> str:
    """
    handle a mismatch we have encountered.
    returns an expression we can tack onto channel (during process / subworkflow call) 
    to fix the mismatch. 
    """
    operations: str = ''
    
    if secondary_secondary_mismatch(srctype, desttype):
        operations += generate_secondary_mismatch_pumbing(srctype, desttype)
    
    elif secondary_single_mismatch(srctype, desttype):
        operations += '.map{ tuple -> tuple[0] }'
    
    if secondary_array_type(desttype):
        operations += '.flatten().toList()'
    
    elif array_type(srctype) and single_type(desttype):
        operations += '.flatten().first()'
    
    elif single_type(srctype) and array_type(desttype):
        operations += '.toList()'
    
    elif secondary_type(srctype) and array_type(desttype):
        operations += '.toList()'
    
    elif secondary_array_type(srctype) and secondary_type(desttype):
        size = get_collate_size(srctype)
        operations += f'.flatten().collate( {size} ).first()'
    
    return operations


# private helpers 
def secondary_secondary_mismatch(srctype: DataType, desttype: DataType) -> bool:
    srctype = nfgen_utils.get_base_type(srctype)
    desttype = nfgen_utils.get_base_type(desttype)
    if secondary_type(srctype) and secondary_type(desttype):
        if srctype.name() != desttype.name():
            return True
    return False

def secondary_single_mismatch(srctype: DataType, desttype: DataType) -> bool:
    # src is secondary, dest is not
    srctype = nfgen_utils.get_base_type(srctype)
    desttype = nfgen_utils.get_base_type(desttype)
    if secondary_type(srctype) and not secondary_type(desttype):
        return True
    # dest is secondary, src is not
    elif secondary_type(desttype) and not secondary_type(srctype):
        return True
    # both secondary or both single
    return False

def generate_secondary_mismatch_pumbing(srctype: DataType, desttype: DataType) -> str:
    """
    1. get the secondary type for srctype & desttype
    2. get the secondary file order for srctype & desttype
    3. iterate through desttype secondary file order, for each find its index in srctype secondary file order
    4. return .map{ tuple -> [tuple[0], tuple[3], tuple[1]] } etc format plumbing
    """
    # 1. get the secondary type for srctype & desttype
    srctype = nfgen_utils.get_base_type(srctype)
    desttype = nfgen_utils.get_base_type(desttype)

    # 2. get the secondary file order for srctype & desttype
    srctype_exts = nfgen_utils.get_extensions(srctype, remove_symbols=True)
    desttype_exts = nfgen_utils.get_extensions(desttype, remove_symbols=True)

    # 3. iterate through desttype secondary file order, for each find its index in srctype secondary file order
    indices: list[int] = []
    for ext in desttype_exts:
        index = srctype_exts.index(ext)
        indices.append(index)

    # 4. return .map{ tuple -> [tuple[0], tuple[3], tuple[1]] } etc format plumbing
    tuple_indices = [f'tuple[{index}]' for index in indices]
    tuple_expr = f"[{', '.join(tuple_indices)}]"
    return f'.map{{ tuple -> {tuple_expr} }}'

