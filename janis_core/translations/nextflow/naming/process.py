

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
    """returns variable names for each file of secondary array type"""
    dtype: DataType = inp.intype  # type: ignore
    exts = utils.get_extensions(dtype, remove_prefix_symbols=True)
    exts = [x.replace('.', '_') for x in exts]
    
    # ['bam', 'bai'] -> ['normal_bam', 'normal_bai'] when other tinput 
    # with same secondary datatype exists
    if duplicate_datatype_exists:
        identifier = to_case(inp.id(), case=settings.translate.nextflow.NF_PROCESS_INPUT_CASE)
        exts = [f'{identifier}_{x}' for x in exts]

    return exts

def secondaries_array(inp: TInput, duplicate_datatype_exists: bool=False) -> str:
    """returns variable name for secondary array input"""
    if duplicate_datatype_exists:
        identifier = to_case(inp.id(), case=settings.translate.nextflow.NF_PROCESS_INPUT_CASE)
        return f'{identifier}_flat'
    else:
        dtype: DataType = inp.intype  # type: ignore
        basetype: DataType = utils.get_base_type(dtype)  # type: ignore
        basetype_name = to_case(basetype.name(), case=settings.translate.nextflow.NF_PROCESS_INPUT_CASE)
        return f'{basetype_name}_flat'

def file_pair(inp: TInput) -> list[str]:
    """returns variable names for each pair for file pair type"""
    basename = to_case(inp.id(), case=settings.translate.nextflow.NF_PROCESS_INPUT_CASE)
    pair1 = f'{basename}1'
    pair2 = f'{basename}2'
    return [pair1, pair2]

def file_pair_array(inp: TInput) -> str:
    """returns variable name for file pair array input"""
    identifier = to_case(inp.id(), case=settings.translate.nextflow.NF_PROCESS_INPUT_CASE)
    return f'{identifier}_flat'


