

from janis_core import (
    ToolInput,
    TInput,
)
from janis_core import settings
from janis_core import translation_utils as utils

from ..casefmt import to_case


def generic(inp: ToolInput | TInput) -> str:
    return to_case(inp.id(), case=settings.translate.nextflow.NF_PROCESS_INPUT_CASE)
    
def secondaries(inp: TInput, duplicate_datatype_exists: bool=False) -> list[str]:
    """returns name of each file for File types with secondaries"""
    # src = sources[inp.id()]
    # srctype: DataType = get_src_type(src)

    # # datatype mismatch! get type info from srctype
    # if srctype.name() != desttype.name():
    #     basetype: File = utils.get_base_type(srctype)  # type: ignore 
    
    # # datatype match. get type info from dest
    # else:
    #     basetype: File = utils.get_base_type(desttype)  # type: ignore 
    
    dtype: DataType = inp.intype  # type: ignore
    # basetype: File = utils.get_base_type(dtype) # type: ignore 
    exts = utils.get_extensions(dtype, remove_prefix_symbols=True)
    exts = [x.replace('.', '_') for x in exts]
    
    # ['bam', 'bai'] -> ['normal_bam', 'normal_bai'] when other tinput 
    # with same secondary datatype exists
    if duplicate_datatype_exists:
        identifier = to_case(inp.id(), case=settings.translate.nextflow.NF_PROCESS_INPUT_CASE)
        exts = [f'{identifier}_{x}' for x in exts]

    return exts

def secondaries_array(inp: TInput, duplicate_datatype_exists: bool=False) -> str:
    if duplicate_datatype_exists:
        identifier = to_case(inp.id(), case=settings.translate.nextflow.NF_PROCESS_INPUT_CASE)
        return f'{identifier}_flat'
    else:
        dtype: DataType = inp.intype  # type: ignore
        basetype: DataType = utils.get_base_type(dtype)  # type: ignore
        basetype_name = to_case(basetype.name(), case=settings.translate.nextflow.NF_PROCESS_INPUT_CASE)
        return f'{basetype_name}_flat'

