


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
from janis_core import translation_utils as utils
from ... import channels
from ... import params



### process inputs

def get_process_inputs(tool: CommandTool | PythonTool, sources: dict[str, Any]) -> set[str]:
    """
    get the tool inputs which will become nextflow process inputs
    """
    if settings.translate.nextflow.MODE == 'workflow':
        return get_process_inputs_workflowmode(tool, sources)
    elif settings.translate.nextflow.MODE == 'tool':  # type: ignore
        return get_process_inputs_toolmode(tool, sources)
    else:
        raise RuntimeError

def get_process_inputs_workflowmode(tool: CommandTool | PythonTool, sources: dict[str, Any]) -> set[str]:
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

def get_process_inputs_toolmode(tool: CommandTool | PythonTool, sources: dict[str, Any]) -> set[str]:
    """all CommandTool inputs. in toolmode, we have no greater scope than the tool itself."""
    return set([x.id() for x in tool.inputs()])


### param inputs

def get_param_inputs(tool: CommandTool | PythonTool, sources: dict[str, Any]) -> set[str]:
    """
    get the tool inputs which will be fed values via params
    """
    if settings.translate.nextflow.MODE == 'workflow':
        return get_param_inputs_workflowmode(tool, sources)
    elif settings.translate.nextflow.MODE == 'tool':  # type: ignore
        return get_param_inputs_toolmode(tool, sources)
    raise RuntimeError('DEV: settings.translate.nextflow.MODE must be either "workflow" or "tool"')

def get_param_inputs_workflowmode(tool: CommandTool | PythonTool, sources: dict[str, Any]) -> set[str]:
    """
    get the inputs which are fed (via step inputs) using a non-File type workflow input
    """
    if settings.translate.nextflow.MINIMAL_PROCESS:
        process_input_ids = get_process_inputs(tool, sources)
        param_input_ids = get_param_process_inputs(sources)
        surviving_ids = param_input_ids - process_input_ids
        return surviving_ids    
    raise NotImplementedError

def get_param_inputs_toolmode(tool: CommandTool | PythonTool, sources: dict[str, Any]) -> set[str]:
    """no param inputs for toolmode. """
    return set()


### internal inputs

def get_internal_inputs(tool: CommandTool | PythonTool, sources: dict[str, Any]) -> set[str]:
    """
    get the tool inputs which will not be exposed to the outside world in the nextflow process
    internal_inputs = all_inputs - process_inputs - param_inputs
    """
    all_inputs: list[TInput] = list(tool.inputs_map().values())

    surviving_ids = get_all_input_ids(all_inputs)
    process_inputs = get_process_inputs(tool, sources)
    param_inputs = get_param_inputs(tool, sources)
    
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
        basetype = utils.get_base_type(inp.intype)
        basetype = utils.ensure_single_type(basetype)
        if isinstance(basetype, (File, Directory, Filename)):
            out.add(inp.id())
    return out

def get_optional_input_ids(tinputs: list[TInput]) -> set[str]:
    """get tool inputs (ids) for tool inputs which are optional"""
    out: set[str] = set()
    for tinp in tinputs:
        basetype = utils.get_base_type(tinp.intype)
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
        node = utils.resolve_node(src)
        if isinstance(node, InputNode):
            if channels.exists(node.uuid):
                out.add(tag)
    return out

def get_param_process_inputs(sources: dict[str, Any]) -> set[str]:
    """get tool inputs (ids) which are being fed a value from a channel"""
    out: set[str] = set()
    for tag, src in sources.items():
        node = utils.resolve_node(src)
        if isinstance(node, InputNode):
            if params.exists(node.uuid):
                out.add(tag)
    return out

def get_connection_process_inputs(sources: dict[str, Any]) -> set[str]:
    """get tool inputs (ids) which are being fed a value from a step connection"""
    out: set[str] = set()
    for tag, src in sources.items():
        node = utils.resolve_node(src)
        if isinstance(node, StepNode):
            out.add(tag)
    return out

def get_scatter_process_inputs(sources: dict[str, Any]) -> set[str]:
    """get tool inputs (ids) which are being scattered on"""
    out: set[str] = set()
    for inname, src in sources.items():
        should_scatter = src.source_map[0].should_scatter
        node = utils.resolve_node(src)
        if should_scatter and isinstance(node, InputNode):
            out.add(inname)
    return out

def get_complex_expression_inputs(sources: dict[str, Any]) -> set[str]:
    """get tool inputs (ids) which are being scattered on"""
    out: set[str] = set()
    for inname, src in sources.items():
        node = utils.resolve_node(src)
        if not isinstance(node, InputNode) and not isinstance(node, StepNode):
            out.add(inname)
    return out
