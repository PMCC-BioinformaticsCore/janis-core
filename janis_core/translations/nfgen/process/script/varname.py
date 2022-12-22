
from typing import Optional, Any
from janis_core import ToolInput, TInput, File

from ... import params
from ... import secondaries
from ... import nfgen_utils


def get_nf_variable_name(inp: ToolInput | TInput, process_inputs: set[str], param_inputs: set[str], sources: dict[str, Any]) -> Optional[str]:
    if inp.id() in process_inputs:
        return process_input_varname(inp, process_inputs, param_inputs)
    elif inp.id() in param_inputs:
        return param_input_varname(inp, sources)
    else:
        raise NotImplementedError

def process_input_varname(inp: ToolInput | TInput, process_inputs: set[str], param_inputs: set[str]) -> Optional[str]:
    # data fed via process input
    name: Optional[str] = None
    if inp.id() in process_inputs:
        # secondary files (name mapped to ext of primary file)
        dtype = inp.input_type if isinstance(inp, ToolInput) else inp.intype
        basetype = nfgen_utils.get_base_type(dtype)
        if isinstance(basetype, File) and basetype.has_secondary_files():
            names = secondaries.get_names(basetype)
            name = names[0]
        # everything else
        else:
            name = inp.id()
    return name
    
def param_input_varname(inp: ToolInput | TInput, sources: dict[str, Any]) -> Optional[str]: 
    # data fed via global param
    src = sources[inp.id()]
    sel = src.source_map[0].source
    param = params.get(sel.input_node.uuid)
    return f'params.{param.name}'