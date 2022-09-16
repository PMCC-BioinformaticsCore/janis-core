

from janis_core import File
from janis_core import Array
from janis_core import DataType
from janis_core import ToolInput
from janis_core.workflow.workflow import InputNode



def get_base_type(task_input: ToolInput | InputNode) -> DataType:
    datatype = task_input.input_type if isinstance(task_input, ToolInput) else task_input.datatype
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

def get_references(task_input: ToolInput | InputNode) -> list[str]:
    """returns the tags of each entity referencing task_input"""
    # TODO
    return [task_input.id()]
    #raise NotImplementedError

def is_nullable(task_input: ToolInput | InputNode) -> bool:
    raise NotImplementedError
