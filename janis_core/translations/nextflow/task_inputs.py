

from typing import Optional
from dataclasses import dataclass, field

from janis_core import (
    ToolInput,
    TInput,
    Tool,
)
from enum import Enum, auto


class TaskInputType(Enum): 
    TASK_INPUT  = auto()
    PARAM       = auto()
    STATIC      = auto()
    IGNORED     = auto()
    LOCAL       = auto()


@dataclass 
class TaskInput:
    tinput_id: str
    value: Optional[str | list[str]]
    ti_type: TaskInputType


@dataclass
class TaskInputRegister:
    
    data_structure: dict[str, dict[str, TaskInput]] = field(default_factory=dict)
    
    """ 
    for each Workflow / CommandTool / PythonTool, stores a map of tinput_id: data_source.
    the data_source is the initial value which drives the ToolInput inside a task. 
    
    data_structure description:
    {
        Scope: {
            'my_file': str,
            'my_string': None,
            'my_secondary': list[str],
            'my_secondary_arr': str,
        }
    }
    """
    def add(
        self, 
        tool_id: str,
        dstype_str: str,
        tinput_id: str, 
        value: Optional[str | list[str]],
    ) -> None:
        """
        adds some data to our data_structure. 
        scope: tool scope
        
        """
        # if scope not yet in data structure, add
        if tool_id not in self.data_structure:
            self.data_structure[tool_id] = {}

        # cast the data source subtype to enum value
        subtype_map = {
            'task_input': TaskInputType.TASK_INPUT,
            'param': TaskInputType.PARAM,
            'static': TaskInputType.STATIC,
            'ignored': TaskInputType.IGNORED,
            'local': TaskInputType.LOCAL,
        }

        # create new DataSource & add to data_structure
        dstype = subtype_map[dstype_str]
        ds = TaskInput(tinput_id, value, dstype)
        self.data_structure[tool_id][tinput_id] = ds

    def get(self, tool_id: str, tinput_id: str) -> TaskInput:
        return self.data_structure[tool_id][tinput_id]
    
    def getall(self, tool_id: str) -> list[TaskInput]:
        return list(self.data_structure[tool_id].values())
    
    def to_string(self) -> str:
        out: str = ''
        for scope_label, name_map in self.data_structure.items():
            out += f'\n{scope_label}\n'
            for tinput_id, ds in name_map.items():
                out += f'{tinput_id}: {ds.value}\n'
        return out
    

ti_register = TaskInputRegister()

def exists(tool_id: str, inp: ToolInput | TInput) -> bool:
    # checks the tool id has an entry for this tool input
    if tool_id not in ti_register.data_structure:
        return False
    if inp.id() not in ti_register.data_structure[tool_id]:
        return False
    return True

def existsall(tool: Tool) -> bool:
    # checks the tool id has an entry 
    # checks each tinput id has an associated value in the entry
    if tool.id() not in ti_register.data_structure:
        return False
    for tinput in tool.tool_inputs():
        if tinput.id() not in ti_register.data_structure[tool.id()]:
            return False
    return True

def get(tool_id: str, inp: ToolInput | TInput) -> TaskInput:
    return ti_register.get(tool_id, inp.id())

def getall(tool_id: str) -> list[TaskInput]:
    return ti_register.getall(tool_id)

def update(tool_id: str, dstype_str: str, tinput_id: str, value: Optional[str | list[str]]):
    ti_register.add(tool_id, dstype_str, tinput_id, value)

def task_inputs(tool_id: str) -> set[str]:
    all_inputs = ti_register.getall(tool_id)
    return set([x.tinput_id for x in all_inputs if x.ti_type == TaskInputType.TASK_INPUT])

def param_inputs(tool_id: str) -> set[str]:
    all_inputs = ti_register.getall(tool_id)
    return set([x.tinput_id for x in all_inputs if x.ti_type == TaskInputType.PARAM])

def static_inputs(tool_id: str) -> set[str]:
    all_inputs = ti_register.getall(tool_id)
    return set([x.tinput_id for x in all_inputs if x.ti_type == TaskInputType.STATIC])

def ignored_inputs(tool_id: str) -> set[str]:
    all_inputs = ti_register.getall(tool_id)
    return set([x.tinput_id for x in all_inputs if x.ti_type == TaskInputType.IGNORED])

def local_inputs(tool_id: str) -> set[str]:
    all_inputs = ti_register.getall(tool_id)
    return set([x.tinput_id for x in all_inputs if x.ti_type == TaskInputType.LOCAL])

def clear() -> None:
    ti_register.data_structure = {}

