

import ast
import regex as re
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
    Workflow,
    CommandTool,
    PythonTool
)

from collections import defaultdict
from janis_core.workflow.workflow import InputNode

from janis_core.graph.steptaginput import StepTagInput
from janis_core.operators.operator import IndexOperator
from janis_core.operators.standard import FirstOperator
from janis_core.operators.selectors import InputNodeSelector, StepOutputSelector, AliasSelector



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



### misc helper methods
# def is_channel_input() -> bool:
#     pass


def get_construct_name(tool: CommandTool | PythonTool | Workflow, scope: list[str]) -> str:
    construct_type = ''
    if isinstance(tool, CommandTool) or isinstance(tool, PythonTool):
        construct_type = 'process'
    elif isinstance(tool, Workflow) and not scope:
        construct_type = 'main_workflow'
    elif isinstance(tool, Workflow):
        construct_type = 'sub_workflow'
    else:
        raise NotImplementedError
    return construct_type

def resolve_node(node: Any) -> Any:
    if isinstance(node, StepTagInput):
        return resolve_node(node.source_map[0].source)
    elif isinstance(node, InputNodeSelector):
        return node.input_node
    elif isinstance(node, StepOutputSelector):
        return node.node
    elif isinstance(node, IndexOperator):
        return resolve_node(node.args[0])
    elif isinstance(node, FirstOperator):
        return resolve_node(node.args[0])
    elif isinstance(node, AliasSelector):
        return resolve_node(node.inner_selector)
    elif isinstance(node, list):
        # TODO this is first selector?
        if len(node) > 0:
            return node[0]
        return None
    else:
        raise NotImplementedError

def get_connections(inp: InputNode, wf: Workflow) -> dict[str, list[str]]:
    connected: dict[str, list[str]] = defaultdict(list)
    for step in wf.step_nodes.values():
        for tinp_id, src in step.sources.items():
            sel = src.source_map[0].source
            if isinstance(sel, InputNodeSelector) and sel.input_node.id() == inp.id():
                connected[step.id()].append(tinp_id)
    return connected

def items_with_id(the_list: list[Any], ids: set[str]) -> list[Any]:
    return [x for x in the_list if x.id() in ids]

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

def get_base_type_task_input(task_input: ToolInput | InputNode | TInput) -> DataType:
    match task_input:
        case ToolInput():
            dtype = task_input.input_type
        case InputNode():
            dtype = task_input.datatype
        case TInput():
            dtype = task_input.intype
        case _:
            raise NotImplementedError
    return get_base_type(dtype)

def is_path(task_input: ToolInput | InputNode) -> bool:
    datatype = get_base_type_task_input(task_input)
    if isinstance(datatype, File):
        return True
    return False

def is_file_pair(task_input: ToolInput | InputNode) -> bool:
    datatype = get_base_type_task_input(task_input)
    if isinstance(datatype, File):
        if datatype.has_secondary_files():
            if len(datatype.secondary_files()) == 1:
                return True
            if len(datatype.secondary_files()) > 1:
                raise NotImplementedError(f'{task_input.id()} has multiple secondaries!')
    return False

def is_nullable(task_input: ToolInput | InputNode) -> bool:
    raise NotImplementedError

def is_simple_path(text: str) -> bool:
    PATH = r'[\w./]+'
    text = text.strip('\'"')
    matches = re.finditer(PATH, text)
    for m in matches:
        if m.span()[1] - m.span()[0] == len(text):
            return True
    return False

def to_groovy(val: Any, dtype: Optional[DataType]=None) -> Any:
    # must work with str version. 
    dtype = get_base_type(dtype)
    val = str(val)
    
    # # secondary files
    # if val != 'None' and isinstance(dtype, File) and dtype.has_secondary_files():
    #     primary_file = wrap(val)
    #     secondary_files: list[str] = []
    #     for suffix in dtype.secondary_files():
    #         sec_file = apply_secondary_file_format_to_filename(primary_file, suffix)
    #         sec_file = wrap(sec_file)
    #         secondary_files.append(sec_file)
    #     # Note: we want primary file to always be the first item in the array
    #     val = [primary_file] + secondary_files

    # wrap in quotes if necessary (for actual string values, file paths etc)
    if should_wrap(val, dtype):
        val = wrap(val)

    # remove dollar variable references (unsure if needed)
    if '$' in val:
        val = val.replace('$', '')

    # 'None' -> 'null' etc
    val = cast_keywords(val)
    return val

def should_wrap(val: str, dtype: Optional[DataType]) -> bool:
    # don't quote lists
    try:
        literal_val = ast.literal_eval(val)
        if isinstance(literal_val, list):
            return False
    except Exception as e: 
        pass

    # don't quote None
    if val == 'None':
        return False

    # don't quote outer array, boolean, numeric types
    no_quote_types: list[type[DataType]] = [Boolean, Int, Float]
    if dtype:
        if type(dtype) in no_quote_types:
            return False

    # don't quote nextflow input channel
    if val.startswith('ch_'):
        return False
    
    # don't quote nextflow referenced variable
    if val.startswith('$'):
        return False
    
    # quote everything else
    return True

def wrap(val: Any) -> Any:
    return f"'{val}'"

def cast_keywords(val: str) -> str:
    # this is done in string world - need a better way of handling lists!
    keyword_map: dict[str, str] = {
        'None': 'null',
        'False': 'false',
        'True': 'true',
    }
    for python_val, groovy_val in keyword_map.items():
        if python_val in val:
            val = val.replace(python_val, groovy_val)
    return val

def get_base_type(dtype: Optional[DataType]) -> Optional[DataType]:
    if dtype is None:
        return dtype
    while isinstance(dtype, Array):
        dtype = dtype.subtype()
    return dtype





