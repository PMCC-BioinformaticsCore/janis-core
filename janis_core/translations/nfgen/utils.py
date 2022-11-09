

from typing import Any, Optional
from janis_core.graph.node import Node
from janis_core.graph.steptaginput import StepTagInput
from janis_core.types import (
    File,
    Boolean,
    Int,
    Float,
    Array
)
from janis_core import (
    DataType,
    ToolInput,
    TInput,
    CommandTool,
    Workflow
)

from janis_core.workflow.workflow import InputNode, StepNode
from janis_core.utils.secondary import apply_secondary_file_format_to_filename
from janis_core.translations.nfgen import ordering
from janis_core.translations.nfgen import settings
import regex as re




"""
FOR THIS SECTION

In the current approach, only Files are supplied using channels. 
Process inputs can therefore only be File types. 
All other inputs are supplied a value using params. 

Additionally, when doing MINIMAL_PROCESS tool translation:
    - inputs with a default are autofilled
    - inputs which are optional are ignored

If doing workflow translation, we get extra info about process inputs
using step inputs. 

For example in MINIMAL_PROCESS mode, some optional files may be being passed 
values in the workflow. these would usually be ignored due to their optionality, but here 
shoud should be included.
"""



### workflow inputs

def get_channel_input_ids(wf: Workflow) -> set[str]:
    """
    Get the (assumed) true workflow inputs. 
    Assume that a workflow input is an InputNode which:
        - has the 'File' datatype
        - is referenced in a step input
    Everything else are static step inputs, or non-exposed tool inputs. 
    """

    # wf inputs with reference in step 
    referenced_inputs = get_referenced_wf_inputs(wf)

    # wf inputs with file type
    file_inputs = get_file_wf_inputs(wf)

    channel_inputs: list[InputNode] = []
    for name, inp in wf.input_nodes.items():
        if name in referenced_inputs and name in file_inputs:
            channel_inputs.append(inp)
    
    # final ordering
    channel_inputs = ordering.workflow_inputs(channel_inputs)
    return {x.id() for x in channel_inputs}

def get_referenced_wf_inputs(wf: Workflow) -> set[str]:
    out: set[str] = set()
    for step in wf.step_nodes.values():
        for src in step.sources.values():
            node = resolve_node(src)
            if isinstance(node, InputNode):
                out.add(node.identifier)
    return out
    
def get_file_wf_inputs(wf: Workflow) -> set[str]:
    out: set[str] = set()
    for name, inp in wf.input_nodes.items():
        dtype = get_base_type(inp)
        if isinstance(dtype, File):
            out.add(name)
    return out



### process inputs

def get_process_input_ids(tool: CommandTool, sources: dict[str, Any]) -> set[str]:
    """
    determine the tool inputs which should remnain as process inputs
    """
    if settings.MODE == 'workflow':
        return get_process_inputs_workflowmode(tool, sources)
    elif settings.MODE == 'tool':  # type: ignore
        return get_process_inputs_toolmode(tool, sources)
    else:
        raise RuntimeError

def get_process_inputs_workflowmode(tool: CommandTool, sources: dict[str, Any]) -> set[str]:
    """
    inputs which are fed (via step inputs) using a file type workflow input
    inputs which are fed (via step inputs) using a connection
    """    
    file_wfinp_ids = get_file_wfinp_connected_input_ids(sources)
    step_conn_ids = get_step_connected_input_ids(sources)
    surviving_ids = file_wfinp_ids | step_conn_ids
    return surviving_ids

def get_process_inputs_toolmode(tool: CommandTool, sources: dict[str, Any]) -> set[str]:
    """
    inputs which are file types
    non-file types usually fed values from params instead.
    
    if MINIMAL_PROCESS:
        - remove inputs which are optional
        - remove inputs with defaults
    """
    all_inputs = tool.inputs()
    
    surviving_ids = get_all_input_ids(all_inputs)
    file_ids = get_file_input_ids(all_inputs)
    optional_ids = get_optional_input_ids(all_inputs)
    default_ids = get_default_input_ids(all_inputs)
    
    if settings.MINIMAL_PROCESS:
        surviving_ids = surviving_ids & file_ids
        surviving_ids = surviving_ids - optional_ids
        surviving_ids = surviving_ids - default_ids
    else:
        surviving_ids = surviving_ids & file_ids

    return surviving_ids


### param inputs

def get_param_input_ids(tool: CommandTool, sources: dict[str, Any]) -> set[str]:
    """
    determine the tool inputs which should be fed a value via params
    """
    if settings.MODE == 'workflow':
        return get_param_inputs_workflowmode(tool, sources)
    elif settings.MODE == 'tool':  # type: ignore
        return get_param_inputs_toolmode(tool, sources)
    else:
        raise RuntimeError

def get_param_inputs_workflowmode(tool: CommandTool, sources: dict[str, Any]) -> set[str]:
    """
    inputs which are fed (via step inputs) using a non-File type workflow input
    """
    if settings.MINIMAL_PROCESS:
        surviving_ids = get_nonfile_wfinp_connected_input_ids(sources)
    else:
        all_inputs = tool.inputs()
        all_ids = get_all_input_ids(all_inputs)
        process_ids = get_process_input_ids(tool, sources)
        surviving_ids = all_ids - process_ids
    return surviving_ids

def get_param_inputs_toolmode(tool: CommandTool, sources: dict[str, Any]) -> set[str]:
    """
    nonfile types 
    
    if MINIMAL_PROCESS:
        - remove inputs which are optional
        - remove inputs with defaults
    """
    all_inputs = tool.inputs()
    
    surviving_ids = get_all_input_ids(all_inputs)
    file_ids = get_file_input_ids(all_inputs)
    optional_ids = get_optional_input_ids(all_inputs)
    default_ids = get_default_input_ids(all_inputs)
    
    if settings.MINIMAL_PROCESS:
        surviving_ids = surviving_ids - file_ids
        surviving_ids = surviving_ids - optional_ids
        surviving_ids = surviving_ids - default_ids
    else:
        surviving_ids = surviving_ids - file_ids

    return surviving_ids


### internal inputs

def get_internal_input_ids(tool: CommandTool, sources: dict[str, Any]) -> set[str]:
    """
    internal inputs = all inputs - process inputs - param inputs
    """
    all_inputs = tool.inputs()

    surviving_ids = get_all_input_ids(all_inputs)
    process_inputs = get_process_input_ids(tool, sources)
    param_inputs = get_param_input_ids(tool, sources)
    
    surviving_ids = surviving_ids - process_inputs
    surviving_ids = surviving_ids - param_inputs
    return surviving_ids



###  sets helper methods

def get_all_input_ids(tinputs: list[ToolInput]) -> set[str]:
    return {x.id() for x in tinputs}

def get_file_input_ids(tinputs: list[ToolInput]) -> set[str]:
    """get tool inputs (ids) for tool inputs which are File types"""
    return {x.id() for x in tinputs if isinstance(get_base_type(x), File)}

def get_optional_input_ids(tinputs: list[ToolInput]) -> set[str]:
    """get tool inputs (ids) for tool inputs which are optional"""
    return {x.id() for x in tinputs if get_base_type(x).optional == True}

def get_default_input_ids(tinputs: list[ToolInput]) -> set[str]:
    """get tool inputs (ids) for tool inputs with a default value"""
    return {x.id() for x in tinputs if x.default != None}

# def get_non_fed_input_ids(tinputs: list[ToolInput], sources: dict[str, Any]) -> set[str]:
#     return {x.id() for x in tinputs if x.id() not in sources}

def get_file_wfinp_connected_input_ids(sources: dict[str, Any]) -> set[str]:
    """get tool inputs (ids) which are being fed a value from a File type workflow input"""
    out: set[str] = set()
    for inname, src in sources.items():
        node = resolve_node(src)
        if isinstance(node, InputNode):
            if isinstance(get_base_type(node), File):
                out.add(inname)
    return out

def get_nonfile_wfinp_connected_input_ids(sources: dict[str, Any]) -> set[str]:
    """get tool inputs (ids) which are being fed a value from a non-File type workflow input"""
    out: set[str] = set()
    for inname, src in sources.items():
        node = resolve_node(src)
        if isinstance(node, InputNode):
            if not isinstance(get_base_type(node), File):
                out.add(inname)
    return out

def get_step_connected_input_ids(sources: dict[str, Any]) -> set[str]:
    """get tool inputs (ids) which are being fed a value from a step connection"""
    out: set[str] = set()
    for inname, src in sources.items():
        node = resolve_node(src)
        if isinstance(node, StepNode):
            out.add(inname)
    return out



### misc helper methods

def items_with_id(the_list: list[Any], ids: set[str]) -> list[Any]:
    return [x for x in the_list if x.id() in ids]

def resolve_node(source: StepTagInput) -> Optional[Node]:
    # workflow input
    if hasattr(source.source_map[0].source, 'input_node'):
        node = source.source_map[0].source.input_node
    
    # connection
    elif hasattr(source.source_map[0].source, 'node'):
        node = source.source_map[0].source.node
    
    # workflow input / connection -> array type
    elif hasattr(source.source_map[0].source, 'args'):
        upstream = source.source_map[0].source.args[0]
        if hasattr(upstream, 'node'):
            node = upstream.node
        elif hasattr(upstream, 'input_node'):
            node = upstream.input_node
    else:
        node = None
    return node

def get_source_value(source: Any) -> Any:
    node = resolve_node(source)
    if isinstance(node, InputNode):
        return node.default
    return None

def roughly_equivalent(val1: Any, val2: Any) -> bool:
    equivalents = {
        '': ' ',
    }
    map_fwd = equivalents
    map_rev = {v: k for k, v in equivalents.items()}
    if val1 is None and val2 is None:
        return True
    elif str(val1) == str(val2):
        return True
    elif val1 in map_fwd and map_fwd[val1] == val2:
        return True
    elif val1 in map_rev and map_rev[val1] == val2:
        return True
    return False

def get_base_type(task_input: ToolInput | InputNode | TInput) -> DataType:
    match task_input:
        case ToolInput():
            dtype = task_input.input_type
        case InputNode():
            dtype = task_input.datatype
        case _:
            dtype = task_input.intype
    while isinstance(dtype, Array):
        dtype = dtype.subtype()
    return dtype

def is_path(task_input: ToolInput | InputNode) -> bool:
    datatype = get_base_type(task_input)
    if isinstance(datatype, File):
        return True
    return False

def is_file_pair(task_input: ToolInput | InputNode) -> bool:
    datatype = get_base_type(task_input)
    if isinstance(datatype, File):
        if datatype.has_secondary_files():
            if len(datatype.secondary_files()) == 1:
                return True
            if len(datatype.secondary_files()) > 1:
                raise NotImplementedError(f'{task_input.id()} has multiple secondaries!')
    return False

def is_nullable(task_input: ToolInput | InputNode) -> bool:
    raise NotImplementedError

type_keyword_map: dict[str, str] = {
    'None': 'null',
    'False': 'false',
    'True': 'true',
}

def is_simple_path(text: str) -> bool:
    PATH = r'[\w./]+'
    text = text.strip('\'"')
    matches = re.finditer(PATH, text)
    for m in matches:
        if m.span()[1] - m.span()[0] == len(text):
            return True
    return False

def wrap_value(val: Any, inp: Optional[ToolInput | InputNode]):
    """
    val is either the inp.default, or can be an override in the case
    of step inputs. we use the val itself + the inp datatype to know
    how to wrap.
    """
    # get dtype
    dtype = None
    if isinstance(inp, InputNode):
        dtype = inp.datatype
    elif isinstance(inp, ToolInput):
        dtype = inp.input_type
    
    # ensure string
    val = str(val)

    # wrap in quotes if necessary (for actual string values, file paths etc)
    if should_wrap(val, inp):
        val = f'"{val}"'
    
    # secondary files
    if isinstance(dtype, File) and dtype.has_secondary_files():
        primary_file = val
        secondary_files: list[str] = []
        for suffix in dtype.secondary_files():
            sec_file = apply_secondary_file_format_to_filename(primary_file, suffix)
            sec_file = f'"{sec_file}"'
            secondary_files.append(sec_file)
        # Note: we want primary file to always be the first item in the array
        val = [primary_file] + secondary_files

    # cast to correct syntax
    if val in type_keyword_map:
        val = type_keyword_map[val]

    # remove dollar variable references (unsure if needed)
    val = val.replace('$', '')
    return val

def should_wrap(val: str, tool_input: Optional[ToolInput | InputNode]):
    # if its bool or null or channel or list or starts with dollar, don't wrap
    # otherwise, wrap
    no_wrap_types: list[type[DataType]] = [Boolean, Int, Float]
    
    # true, false, numeric
    if tool_input:
        dtype = get_base_type(tool_input)
        if type(dtype) in no_wrap_types:
            return False
        elif dtype.is_array():
            return False
        elif val == 'None':
            return False

    # input channel
    if val.startswith('ch_'):
        return False
    
    # referenced variable
    if val.startswith('$'):
        return False
    
    return True






# def get_exposed_tool_inputs(tool: CommandTool, sources: dict[str, Any]) -> list[ToolInput]:
#     """
#     Tool inputs which rely on outside world
#     Note: tool input must be fed a value from 'sources' (step inputs) to be considered.
    
#     true conditions: 
#     1. source is connection
#     2. source is workflow input
#         - && tool input has file type
#         - && wf input has file type
#         or
#         - && tool input & wf input have same type
#         - && tool input default == wf input default

#     """
#     out: list[ToolInput] = []
#     inputs: list[ToolInput] = tool.inputs()
#     for toolinp in inputs:
#         # must have source (be in step inputs)
#         if toolinp.id() in sources:
#             source = sources[toolinp.id()]
#             node = resolve_node(source)
#             if isinstance(node, StepNode):
#                 out.append(toolinp)
#             elif isinstance(node, InputNode):
#                 toolinp_dtype = get_base_type(toolinp)
#                 wfinp_dtype = get_base_type(node)
#                 if isinstance(toolinp_dtype, File) and isinstance(wfinp_dtype, File):
#                     out.append(toolinp)
#                 elif type(toolinp_dtype) == type(wfinp_dtype):
#                     if not roughly_equivalent(node.default, toolinp.default):
#                         out.append(toolinp)
#     return out
