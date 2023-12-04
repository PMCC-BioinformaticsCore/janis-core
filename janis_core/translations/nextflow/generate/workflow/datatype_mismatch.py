

from janis_core.types import DataType

from janis_core import translation_utils as utils
from janis_core.translation_utils import DTypeType

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
    while dtype.is_array() and dtype.name() == 'Array' and dtype.subtype():
        depth += 1
        dtype = dtype.subtype()
    return depth


# public
def gen_datatype_mismatch_plumbing(srctype: DataType, desttype: DataType, destscatter: bool, is_connection: bool) -> str:
    """
    handle a mismatch we have encountered.
    returns an expression we can tack onto channel (during process / subworkflow call) 
    to fix the mismatch. 
    """
    operations: str = ''
    if is_datatype_mismatch(srctype, desttype, destscatter):
        operations += gen_base_datatype_plumbing(srctype, desttype, is_connection)  # must be first
        operations += gen_array_datatype_plumbing(srctype, desttype, destscatter)
    return operations

def gen_base_datatype_plumbing(srctype: DataType, desttype: DataType, is_connection: bool) -> str:
    """handle base datatype mismatch transformations"""
    base_srctype = utils.get_base_type(srctype)
    base_desttype = utils.get_base_type(desttype)

    src_dtt = utils.get_dtt(base_srctype)
    dest_dtt = utils.get_dtt(base_desttype)

    # secondary -> single
    if src_dtt == DTypeType.SECONDARY and dest_dtt == DTypeType.FILE:
        return generate_secondary_single_plumbing(base_srctype, is_connection)

    # secondary -> secondary
    elif src_dtt == dest_dtt == DTypeType.SECONDARY:
        if base_srctype.name() != base_desttype.name():
            return generate_secondary_secondary_mismatch_pumbing(base_srctype, base_desttype, is_connection)
    
    return ''

def generate_secondary_single_plumbing(srctype: DataType, is_connection: bool) -> str:
    """
    generates plumbing for secondary -> single mismatch
    eg BamBai -> Bam
    """
    if utils.datatype_will_be_channel(srctype) or is_connection:
        return '.map{ tuple -> tuple[0] }'
    elif utils.datatype_will_be_variable(srctype):
        return '[0]'
    else:
        raise RuntimeError(f"Disallowed datatype for secondary -> single plumbing: {srctype}")

def generate_secondary_secondary_mismatch_pumbing(srctype: DataType, desttype: DataType, is_connection: bool) -> str:
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
    for i, ext in enumerate(desttype_exts):
        if ext == 'primary':
            index = 0
        elif ext in srctype_exts:
            index = srctype_exts.index(ext)
        else:
            # print(f"Secondary tuple mapping failed for {ext}. Made a best guess.")
            index = i

        indices.append(index)
    
    if utils.datatype_will_be_channel(srctype) or is_connection:
        # var.map{ tuple -> [tuple[0], tuple[3], tuple[1]] } etc format plumbing
        tuple_indices = [f'tuple[{index}]' for index in indices]
        tuple_expr = f"[{', '.join(tuple_indices)}]"
        return f'.map{{ tuple -> {tuple_expr} }}'
    
    elif utils.datatype_will_be_variable(srctype):
        # var[0, 1, 2] etc format plumbing
        index_list = ', '.join([str(index) for index in indices])
        return f'[{index_list}]'
    
    else:
        raise RuntimeError(f"Disallowed datatypes for secondary -> secondary plumbing: {srctype}, {desttype}")

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


