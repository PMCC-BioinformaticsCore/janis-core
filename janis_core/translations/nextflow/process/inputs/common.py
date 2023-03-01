


from typing import Any

from janis_core.workflow.workflow import InputNode, StepNode
from janis_core import (
    File,
    Filename,
    Directory,
    CommandTool,
    PythonTool,
    TInput,
)

from janis_core import settings
from ... import channels
from ... import params
from ... import nfgen_utils as nfgen_utils



### process inputs

def get_process_inputs(sources: dict[str, Any]) -> set[str]:
    """
    get the tool inputs which will become nextflow process inputs
    """
    if settings.translate.nextflow.MODE == 'workflow':
        return get_process_inputs_workflowmode(sources)
    elif settings.translate.nextflow.MODE == 'tool':  # type: ignore
        return get_process_inputs_toolmode(sources)
    else:
        raise RuntimeError

def get_process_inputs_workflowmode(sources: dict[str, Any]) -> set[str]:
    """
    inputs which are fed (via step inputs) using a file type workflow input
    inputs which are fed (via step inputs) using a connection
    inputs which are involved in scatter
    """
    channel_wfinp_ids = get_channel_process_inputs(sources)
    step_conn_ids = get_connection_process_inputs(sources)
    scatter_wfinp_ids = get_scatter_process_inputs(sources)
    complex_expr_ids = get_complex_expression_inputs(sources)
    surviving_ids = channel_wfinp_ids | step_conn_ids | scatter_wfinp_ids | complex_expr_ids
    return surviving_ids

def get_process_inputs_toolmode(sources: dict[str, Any]) -> set[str]:
    """
    inputs which are file types
    non-file types usually fed values from params instead.
    
    if MINIMAL_PROCESS:
        - remove inputs which are optional
        - remove inputs with defaults
    """
    raise NotImplementedError
    # all_inputs: list[TInput] = list(tool.inputs_map().values())
    
    # surviving_ids = get_all_input_ids(all_inputs)
    # file_ids = get_file_input_ids(all_inputs)
    # optional_ids = get_optional_input_ids(all_inputs)
    # default_ids = get_default_input_ids(all_inputs)
    
    # if settings.MINIMAL_PROCESS:
    #     surviving_ids = surviving_ids & file_ids
    #     surviving_ids = surviving_ids - optional_ids
    #     surviving_ids = surviving_ids - default_ids
    # else:
    #     surviving_ids = surviving_ids & file_ids

    # return surviving_ids




### param inputs

def get_param_inputs(sources: dict[str, Any]) -> set[str]:
    """
    get the tool inputs which will be fed values via params
    """
    if settings.translate.nextflow.MODE == 'workflow':
        return get_param_inputs_workflowmode(sources)
    elif settings.translate.nextflow.MODE == 'tool':  # type: ignore
        return get_param_inputs_toolmode(sources)
    raise RuntimeError('DEV: settings.translate.nextflow.MODE must be either "workflow" or "tool"')

def get_param_inputs_workflowmode(sources: dict[str, Any]) -> set[str]:
    """
    get the inputs which are fed (via step inputs) using a non-File type workflow input
    """
    if settings.translate.nextflow.MINIMAL_PROCESS:
        process_input_ids = get_process_inputs(sources)
        param_input_ids = get_param_process_inputs(sources)
        surviving_ids = param_input_ids - process_input_ids
        return surviving_ids    
    raise NotImplementedError

def get_param_inputs_toolmode(sources: dict[str, Any]) -> set[str]:
    """
    nonfile types 
    
    if MINIMAL_PROCESS:
        - remove inputs which are optional
        - remove inputs with defaults
    """
    raise NotImplementedError
    all_inputs: list[TInput] = list(tool.inputs_map().values())
    
    surviving_ids = get_all_input_ids(all_inputs)
    file_ids = get_file_input_ids(all_inputs)
    optional_ids = get_optional_input_ids(all_inputs)
    default_ids = get_default_input_ids(all_inputs)
    
    if settings.translate.nextflow.MINIMAL_PROCESS:
        surviving_ids = surviving_ids - file_ids
        surviving_ids = surviving_ids - optional_ids
        surviving_ids = surviving_ids - default_ids
    else:
        surviving_ids = surviving_ids - file_ids

    return surviving_ids


### internal inputs

def get_internal_inputs(tool: CommandTool | PythonTool, sources: dict[str, Any]) -> set[str]:
    """
    get the tool inputs which will not be exposed to the outside world in the nextflow process
    internal_inputs = all_inputs - process_inputs - param_inputs
    """
    all_inputs: list[TInput] = list(tool.inputs_map().values())

    surviving_ids = get_all_input_ids(all_inputs)
    process_inputs = get_process_inputs(sources)
    param_inputs = get_param_inputs(sources)
    
    surviving_ids = surviving_ids - process_inputs
    surviving_ids = surviving_ids - param_inputs
    return surviving_ids



###  local helper methods

# based on TInput info
def get_all_input_ids(tinputs: list[TInput]) -> set[str]:
    return {x.id() for x in tinputs}

def get_file_input_ids(tinputs: list[TInput]) -> set[str]:
    """get tool inputs (ids) for tool inputs which are File types"""
    out: set[str] = set()
    for inp in tinputs:
        basetype = nfgen_utils.get_base_type(inp.intype)
        basetype = nfgen_utils.ensure_single_type(basetype)
        if isinstance(basetype, (File, Directory, Filename)):
            out.add(inp.id())
    return out

def get_optional_input_ids(tinputs: list[TInput]) -> set[str]:
    """get tool inputs (ids) for tool inputs which are optional"""
    out: set[str] = set()
    for tinp in tinputs:
        basetype = nfgen_utils.get_base_type(tinp.intype)
        if basetype and basetype.optional:
            out.add(tinp.id())
    return out

def get_default_input_ids(tinputs: list[TInput]) -> set[str]:
    """get tool inputs (ids) for tool inputs with a default value"""
    return {x.id() for x in tinputs if x.default is not None}


# based on step info
def get_channel_process_inputs(sources: dict[str, Any]) -> set[str]:
    """get tool inputs (ids) which are being fed a value from a channel"""
    out: set[str] = set()
    for tag, src in sources.items():
        node = nfgen_utils.resolve_node(src)
        if isinstance(node, InputNode):
            if channels.exists(node.uuid):
                out.add(tag)
    return out

def get_param_process_inputs(sources: dict[str, Any]) -> set[str]:
    """get tool inputs (ids) which are being fed a value from a channel"""
    out: set[str] = set()
    for tag, src in sources.items():
        node = nfgen_utils.resolve_node(src)
        if isinstance(node, InputNode):
            if params.exists(node.uuid):
                out.add(tag)
    return out

def get_connection_process_inputs(sources: dict[str, Any]) -> set[str]:
    """get tool inputs (ids) which are being fed a value from a step connection"""
    out: set[str] = set()
    for tag, src in sources.items():
        node = nfgen_utils.resolve_node(src)
        if isinstance(node, StepNode):
            out.add(tag)
    return out

def get_scatter_process_inputs(sources: dict[str, Any]) -> set[str]:
    """get tool inputs (ids) which are being scattered on"""
    out: set[str] = set()
    for inname, src in sources.items():
        scatter = src.source_map[0].scatter
        node = nfgen_utils.resolve_node(src)
        if scatter and isinstance(node, InputNode):
            out.add(inname)
    return out

def get_complex_expression_inputs(sources: dict[str, Any]) -> set[str]:
    """get tool inputs (ids) which are being scattered on"""
    out: set[str] = set()
    for inname, src in sources.items():
        node = nfgen_utils.resolve_node(src)
        if not isinstance(node, InputNode) and not isinstance(node, StepNode):
            out.add(inname)
    return out
