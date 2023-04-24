

from janis_core.types import DataType

from janis_core import translation_utils as utils

from .common import array_type
from .common import single_type
from .common import secondary_type
from .common import secondary_array_type
from .common import get_collate_size


# public
def is_datatype_mismatch(srctype: DataType, desttype: DataType, destscatter: bool) -> bool:
    if secondary_array_type(desttype):
        # Array(Secondary()) types are always considered datatype mismatch
        # because they get flattened before being fed to a process. 
        return True
    elif is_array_depth_mismatch(srctype, desttype, destscatter):
        return True
    elif is_base_type_mismatch(srctype, desttype):
        return True
    return False

def is_base_type_mismatch(srctype: DataType, desttype: DataType) -> bool:
    base_srctype = utils.get_base_type(srctype)
    base_desttype = utils.get_base_type(desttype)
    if base_srctype.name() != base_desttype.name():
        return True
    return False

def is_array_depth_mismatch(srctype: DataType, desttype: DataType, destscatter: bool) -> bool:
    """
    identify whether the datatypes have array differences.
    eg Array(String()) -> String() 
    """
    srctype_depth = get_array_depth(srctype)
    desttype_depth = get_array_depth(desttype)

    # if srctype array depth equals desttype, no mismatch
    if abs(srctype_depth - desttype_depth) == 0:
        return False
    
    # if srctype array depth is one or more different to desttype, mismatch
    elif abs(srctype_depth - desttype_depth) >= 1:
        return True

    return False

def get_array_depth(dtype: DataType) -> int:
    depth = 0
    while dtype.is_array() and dtype.subtype() and not utils.is_file_pair_type(dtype, recursive=False):
        depth += 1
        dtype = dtype.subtype()
    return depth


# public
def gen_datatype_mismatch_plumbing(srctype: DataType, desttype: DataType, destscatter: bool) -> str:
    """
    handle a mismatch we have encountered.
    returns an expression we can tack onto channel (during process / subworkflow call) 
    to fix the mismatch. 
    """
    operations: str = ''
    if is_datatype_mismatch(srctype, desttype, destscatter):
        operations += gen_base_datatype_plumbing(srctype, desttype)  # must be first
        operations += gen_array_datatype_plumbing(srctype, desttype, destscatter)
    return operations

def gen_base_datatype_plumbing(srctype: DataType, desttype: DataType) -> str:
    """handle base datatype mismatch transformations"""
    base_srctype = utils.get_base_type(srctype)
    base_desttype = utils.get_base_type(desttype)
    
    if secondary_secondary_mismatch(base_srctype, base_desttype):
        return generate_secondary_mismatch_pumbing(srctype, desttype)
    
    elif secondary_single_mismatch(base_srctype, base_desttype):
        return '.map{ tuple -> tuple[0] }'
    
    return ''

def gen_array_datatype_plumbing(srctype: DataType, desttype: DataType, destscatter: bool) -> str:
    """handle array depth mismatch transformations"""
    srctype_depth = get_array_depth(srctype)
    desttype_depth = get_array_depth(desttype)

    if secondary_array_type(desttype):
        return '.flatten().toList()'

    # if srctype array depth is 1 and desttype is 0, we use flatten()
    if srctype_depth == 1 and desttype_depth == 0:
        # ([bam, bai]) -> (bam)
        if secondary_array_type(srctype) and single_type(desttype):
            return f'.flatten().first()'
        
        # ([[bam, bai], [bam, bai]]) -> ([bam, bai], [bam, bai])
        elif secondary_array_type(srctype) and secondary_type(desttype) and destscatter:
            size = get_collate_size(desttype)
            return f'.flatten().collate( {size} )'
        
        # ([[bam, bai], [bam, bai]]) -> ([bam, bai])
        elif secondary_array_type(srctype) and secondary_type(desttype) and not destscatter:
            size = get_collate_size(desttype)
            return f'.flatten().collate( {size} ).first()'
        
        # ([bam, bam]) -> (bam, bam)
        elif array_type(srctype) and single_type(desttype) and destscatter:
            return '.flatten()'
        
        # ([bam]) -> (bam)
        elif array_type(srctype) and single_type(desttype) and not destscatter:
            return '.flatten().first()'
    
    # if srctype array depth is 0 and desttype is 1, we use toList() 
    elif srctype_depth == 0 and desttype_depth == 1:
        return '.toList()'

    elif srctype_depth >= 2 or desttype_depth >= 2:
        raise NotImplementedError('2 or more levels of array nesting')

    return ''



# private helpers below ----
def secondary_secondary_mismatch(srctype: DataType, desttype: DataType) -> bool:
    if secondary_type(srctype) and secondary_type(desttype):
        if srctype.name() != desttype.name():
            return True
    return False

def secondary_single_mismatch(srctype: DataType, desttype: DataType) -> bool:
    # src is secondary, dest is not
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
    srctype = utils.get_base_type(srctype)
    desttype = utils.get_base_type(desttype)

    # 2. get the secondary file order for srctype & desttype
    srctype_exts = utils.get_extensions(srctype, remove_prefix_symbols=True)
    desttype_exts = utils.get_extensions(desttype, remove_prefix_symbols=True)

    # 3. iterate through desttype secondary file order, for each find its index in srctype secondary file order
    indices: list[int] = []
    for ext in desttype_exts:
        index = srctype_exts.index(ext)
        indices.append(index)

    # 4. return .map{ tuple -> [tuple[0], tuple[3], tuple[1]] } etc format plumbing
    tuple_indices = [f'tuple[{index}]' for index in indices]
    tuple_expr = f"[{', '.join(tuple_indices)}]"
    return f'.map{{ tuple -> {tuple_expr} }}'

