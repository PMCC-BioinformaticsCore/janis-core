

from typing import Any, Optional
from janis_core import File
from janis_core import Boolean
from janis_core import Int
from janis_core import Float
from janis_core import Array
from janis_core import DataType
from janis_core import ToolInput
from janis_core import TInput
from janis_core.workflow.workflow import InputNode



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

def wrap_value(val: Any, tool_input: Optional[ToolInput | InputNode]):
    # ensure string
    val = str(val)    
    
    # wrap in quotes if necessary (for actual string values, file paths etc)
    if should_wrap(val, tool_input):
        val = f"'{val}'"
    
    # cast to correct syntax
    if val in type_keyword_map:
        val = type_keyword_map[val]

    # remove dollar variable references (unsure if needed)
    val = val.replace('$', '')
    return val

def should_wrap(val: str, tool_input: Optional[ToolInput | InputNode]):
    # if its bool or null or channel or starts with dollar, don't wrap
    # otherwise, wrap
    no_wrap_types: list[type[DataType]] = [Boolean, Int, Float]
    
    # true, false, numeric
    if tool_input:
        dtype = get_base_type(tool_input)
        if type(dtype) in no_wrap_types:
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
