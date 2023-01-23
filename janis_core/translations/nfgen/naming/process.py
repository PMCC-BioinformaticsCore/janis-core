

from typing import Optional, Any

from janis_core.types import File, Stdout
from janis_core import (
    DataType,
    ToolInput,
    TInput,
)

from ..casefmt import to_case
from ..plumbing import trace
from .. import nfgen_utils
from .. import params
from .. import settings


### PROCESS SPECIFIC

# this should be kept in ProcessVariableManager or something
def get_varname_toolinput(inp: ToolInput | TInput, process_inputs: set[str], param_inputs: set[str], sources: dict[str, Any]) -> Optional[str]:
    if inp.id() in process_inputs:
        return _get_varname_process_input(inp, sources, process_inputs, param_inputs)
    elif inp.id() in param_inputs:
        return _get_varname_param_input(inp, sources)
    else:
        return None
        # internal inputs in the future
        raise NotImplementedError

def _get_varname_param_input(inp: ToolInput | TInput, sources: dict[str, Any]) -> Optional[str]: 
    # data fed via global param
    src = sources[inp.id()]
    sel = src.source_map[0].source
    param = params.get(sel.input_node.uuid)
    return f'params.{param.name}'

def _get_varname_process_input(inp: ToolInput | TInput, sources: dict[str, Any], process_inputs: set[str], param_inputs: set[str]) -> Optional[str]:
    # data fed via process input
    name: Optional[str] = None
    
    if inp.id() in process_inputs:
        dtype: DataType = inp.input_type if isinstance(inp, ToolInput) else inp.intype  # type: ignore
        if nfgen_utils.is_array_secondary_type(dtype):
            name = process_input_secondaries_array(inp)
        elif nfgen_utils.is_secondary_type(dtype):
            name = process_input_secondaries(inp, sources)[0]  # first extension
        else:
            name = process_input_generic(inp)  
    
    return name

def process_input_secondaries_array(inp: ToolInput | TInput) -> str:
    dtype: DataType = inp.input_type if isinstance(inp, ToolInput) else inp.intype  # type: ignore
    basetype: DataType = nfgen_utils.get_base_type(dtype)  # type: ignore
    basetype_name = to_case(basetype.name(), case=settings.NF_PROCESS_INPUT_CASE)
    return f'{basetype_name}_array_flat'

def process_input_secondaries_array_primary_files(inp: ToolInput | TInput) -> str:
    """
    example: ToolInput = Array(BamBai)
    process input name: indexed_bam_array_flat
    primary files name: bams
    """
    dtype: DataType = inp.input_type if isinstance(inp, ToolInput) else inp.intype  # type: ignore
    basetype: File = nfgen_utils.get_base_type(dtype)  # type: ignore 
    exts = nfgen_utils.get_extensions(basetype, remove_symbols=True)
    primary_ext = exts[0]
    primary_name = f'{primary_ext}s' # primary extension -> make plural
    return primary_name

def process_input_secondaries(inp: ToolInput | TInput, sources: dict[str, Any]) -> list[str]:
    """returns name of each file for File types with secondaries"""

    src = sources[inp.id()]
    srctype: DataType = get_src_type(src)
    desttype: DataType = inp.input_type if isinstance(inp, ToolInput) else inp.intype  # type: ignore

    # datatype mismatch! get type info from srctype
    if srctype.name() != desttype.name():
        basetype: File = nfgen_utils.get_base_type(srctype)  # type: ignore 
    
    # datatype match. get type info from dest
    else:
        basetype: File = nfgen_utils.get_base_type(desttype)  # type: ignore 
    
    exts = nfgen_utils.get_extensions(basetype, remove_symbols=True)
    return exts
    
def get_src_type(src: Any) -> DataType:
    # the srctype corresponds to either a workflow input, or step output.
    # scattering doesn't matter. 
    dtype = trace.trace_source_datatype(src)
    if isinstance(dtype, Stdout):
        return dtype.subtype
    else:
        return dtype

def process_input_generic(inp: ToolInput | TInput) -> str:
    return to_case(inp.id(), case=settings.NF_PROCESS_INPUT_CASE)