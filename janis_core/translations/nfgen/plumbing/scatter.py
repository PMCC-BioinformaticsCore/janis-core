


from janis_core.types import DataType

from .. import nfgen_utils
from .. import settings

from .common import array_type
from .common import single_type
from .common import secondary_type
from .common import secondary_array_type

INDENT = settings.NF_INDENT


def is_scatter_relationship(srcscatter: bool, destscatter: bool) -> bool:
    if srcscatter or destscatter:
        return True
    return False

def handle_scatter_relationship(srcscatter: bool, destscatter: bool, srctype: DataType, desttype: DataType) -> str:
    # non-secondaries
    if _scatter_to_scatter(srcscatter, destscatter, srctype, desttype):
        return '.toList().flatten()'
    
    elif _scatter_to_array(srcscatter, destscatter, srctype, desttype):
        return '.toList()'
    
    elif _array_to_scatter(srcscatter, destscatter, srctype, desttype):
        return '.flatten()'
    
    # secondaries
    elif _scatter_secondary_array_to_scatter_secondary(srcscatter, destscatter, srctype, desttype):
        basetype = nfgen_utils.get_base_type(srctype)
        exts = nfgen_utils.get_extensions(basetype)
        size = len(exts)
        return f'.toList().flatten().collate( {size} )'
    
    elif _scatter_secondary_array_to_secondary_array(srcscatter, destscatter, srctype, desttype):
        return '.flatten().toList()'
    
    elif _secondary_array_to_scatter_secondary_array(srcscatter, destscatter, srctype, desttype):
        basetype = nfgen_utils.get_base_type(srctype)
        exts = nfgen_utils.get_extensions(basetype)
        size = len(exts)
        return f'.flatten().collate( {size} )'
    
    else:
        raise RuntimeError('DEV: there should be a scatter relationship here, but apparently not?')


# non-secondaries
def _scatter_to_scatter(srcscatter: bool, destscatter: bool, srctype: DataType, desttype: DataType) -> bool:
    if srcscatter and destscatter:
        if single_type(srctype) and single_type(desttype): 
            return True
    return False

def _scatter_to_array(srcscatter: bool, destscatter: bool, srctype: DataType, desttype: DataType) -> bool:
    # stp1 scattered, stp2 not scattered
    if srcscatter and not destscatter:
        if single_type(srctype) and array_type(desttype):
            return True
    return False

def _array_to_scatter(srcscatter: bool, destscatter: bool, srctype: DataType, desttype: DataType) -> bool:
    if not srcscatter and destscatter:
        if array_type(srctype) and single_type(desttype):
            return True
    return False


# secondaries
def _scatter_secondary_array_to_scatter_secondary(srcscatter: bool, destscatter: bool, srctype: DataType, desttype: DataType) -> bool:
    if srcscatter and destscatter:
        if secondary_type(srctype) and secondary_type(desttype): 
            return True
    return False

def _scatter_secondary_array_to_secondary_array(srcscatter: bool, destscatter: bool, srctype: DataType, desttype: DataType) -> bool:
    # stp1 scattered, stp2 not scattered
    if srcscatter and not destscatter:
        if secondary_type(srctype) and secondary_array_type(desttype):
            return True
    return False

def _secondary_array_to_scatter_secondary_array(srcscatter: bool, destscatter: bool, srctype: DataType, desttype: DataType) -> bool:
    # stp1 not scattered, stp2 scattered
    if not srcscatter and destscatter:
        if secondary_array_type(srctype) and secondary_type(desttype):
            return True
    return False




# SCATTER: CARTESIAN CROSS

CROSS_CHANNEL_NAME = 'ch_cartesian_cross'

def cartesian_cross_subname(channel: str) -> str:
    # eg: ch_bams -> ch_cartesian_cross.bams
    ch_subname = channel.replace('ch_', '')
    return f'{CROSS_CHANNEL_NAME}.{ch_subname}'   

def cartesian_cross_operation(channels: list[str]) -> str:
    """
    * just an example, you wouldnt actually do this cuz it mixes up the bams and bais

    - janis - 
    scatter=ScatterDescription(
        ["bams", "bais"],
        method=ScatterMethod.cross
    ),

    - nextflow - 
    ch_bams
    .combine(ch_bais)
    .multiMap { it ->
        bams: it[0]
        bais: it[1]
    }
    .set { ch_cartesian_cross }
    """
    # .flatten()?????

    lines: list[str] = []
    lines.append(channels[0])                   # ch_bams
    
    # cartesian cross
    for ch_name in channels[1:]:
        lines.append(f'.combine({ch_name})')    # .combine(ch_bais)
    
    # multiMap separation
    lines.append('.multiMap { it ->')           # .multiMap { it ->
    for i, ch_name in enumerate(channels):
        ch_subname = ch_name.replace('ch_', '')     
        lines.append(f'{INDENT}{ch_subname}: it[{i}]')     # bams: it[0] \n bais: it[1]
    lines.append('}')

    # set new channel names
    lines.append(f'.set {{ {CROSS_CHANNEL_NAME} }}')
    text = '\n'.join(lines)
    return text






