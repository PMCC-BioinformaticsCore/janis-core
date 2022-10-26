

from typing import Any, Optional
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
import regex as re


def get_workflow_inputs(wf: Workflow) -> list[InputNode]:
    """
    Get the (assumed) true workflow inputs. 
    Assume that a workflow input is an InputNode which:
        - has the 'File' datatype
        - is referenced in a step input
    Everything else are static step inputs, or non-exposed tool inputs. 
    """
    # wf inputs with reference in step 
    referenced_inputs: set[str] = set()
    for step in wf.step_nodes.values():
        for src in step.sources.values():
            if hasattr(src.source_map[0], 'source'):
                if hasattr(src.source_map[0].source, 'input_node'):
                    referenced_inputs.add(src.source_map[0].source.input_node.identifier)
    
    # get names of true wf inputs using sets
    all_wf_inputs = wf.input_nodes
    file_inputs = set([k for k, v in all_wf_inputs.items() if isinstance(v.datatype, File)])

    # get true_wf_inputs using input names
    out: list[InputNode] = []
    for name, wfinp in all_wf_inputs.items():
        if name in file_inputs and name in referenced_inputs:
            out.append(wfinp)

    # final ordering
    wfinps = ordering.workflow_inputs(out)
    return wfinps

def get_exposed_tool_inputs(tool: CommandTool) -> list[ToolInput]:
    inputs = [x for x in tool.inputs() if x.id() in tool.connections]
    inputs = ordering.tool_inputs(inputs)
    return inputs

def get_channel_tool_inputs(tool: CommandTool, values: dict[str, Any]) -> list[ToolInput]:
    """
    start with all exposed tool inputs
    keep Files
    keep those with wfinput or connection sources
    """
    all_inputs = get_exposed_tool_inputs(tool)
    file_inputs = [x for x in all_inputs if isinstance(x.input_type, File)]
    var_inputs: set[str] = set()
    for inname, src in values.items():
        if hasattr(src, 'source_map'):
            if hasattr(src.source_map[0], 'source'):
                if hasattr(src.source_map[0].source, 'input_node'):
                    var_inputs.add(inname)
                if hasattr(src.source_map[0].source, 'node'):
                    var_inputs.add(inname)
    
    channel_inputs = [x for x in file_inputs if x.id() in var_inputs]
    return channel_inputs

def get_param_tool_inputs(tool: CommandTool, values: dict[str, Any]) -> list[ToolInput]:
    """
    Its better to find NOT channel_inputs, than it is to have custom logic to 
    define param_inputs. This way we are sure each input is either in 
    channel_inputs or param_inputs.
    Logic:
        > all_inputs = channel_inputs + param_inputs
        > param_inputs = all_inputs - channel_inputs
    """
    all_inputs = get_exposed_tool_inputs(tool)
    channel_inputs = get_channel_tool_inputs(tool, values)
    channel_input_names = set([x.id() for x in channel_inputs])
    param_inputs = [x for x in all_inputs if x.id() not in channel_input_names]
    return param_inputs

def get_base_type(task_input: ToolInput | InputNode | TInput) -> DataType:
    match task_input:
        case ToolInput():
            datatype = task_input.input_type
        case InputNode():
            datatype = task_input.datatype
        case _:
            datatype = task_input.intype
    while isinstance(datatype, Array):
        datatype = datatype.subtype()
    return datatype

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

def get_input_references(inp: InputNode) -> list[str]:
    """returns the tags of each entity referencing task_input"""
    # TODO!
    return [inp.id()]
    #raise NotImplementedError

def is_nullable(task_input: ToolInput | InputNode) -> bool:
    raise NotImplementedError


type_keyword_map: dict[str, str] = {
    'None': 'null',
    'False': 'false',
    'True': 'true',
    "''": 'null',
    '': 'null',
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
        val = f"'{val}'"
    
    # secondary files
    if isinstance(dtype, File) and dtype.has_secondary_files():
        primary_file = val
        secondary_files: list[str] = []
        for suffix in dtype.secondary_files():
            sec_file = apply_secondary_file_format_to_filename(primary_file, suffix)
            sec_file = f"'{sec_file}'"
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
