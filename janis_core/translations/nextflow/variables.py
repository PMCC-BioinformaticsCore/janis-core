
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Optional, Any
from enum import Enum, auto
from dataclasses import dataclass, field

from janis_core import Tool, WorkflowBuilder, TInput

from . import task_inputs
from .task_inputs import TaskInputType


def init_variable_manager_for_task(tool: Tool) -> VariableManager:
    """
    initialise a variable manager for a tool, covering all tool inputs. 
    differs from task_inputs because task_inputs are based on the task 'class', whereas the variable manager is specific to an 'instance' of the task 
    """
    if isinstance(tool, WorkflowBuilder):
        initialiser = WorkflowVariableManagerInitialiser(tool)
    else:
        initialiser = ProcessVariableManagerInitialiser(tool)
    return initialiser.initialise()


class VariableManagerInitialiser(ABC):
    types_map = {
        TaskInputType.TASK_INPUT: 'task_input',
        TaskInputType.PARAM: 'param',
        TaskInputType.STATIC: 'static',
        TaskInputType.IGNORED: 'ignored',
        TaskInputType.LOCAL: 'local',
    }

    def __init__(self, tool: Tool) -> None:
        self.tool = tool

    @abstractmethod
    def initialise(self) -> VariableManager:
        ...


class ProcessVariableManagerInitialiser(VariableManagerInitialiser):
    
    def initialise(self) -> VariableManager:
        vmanager = VariableManager()
        
        for tinput in self.tool.tool_inputs():
            task_input = task_inputs.get(self.tool.id(), tinput)
            vtype_str = self.types_map[task_input.ti_type]
            vmanager.update(task_input.tinput_id, vtype_str=vtype_str, value=task_input.value)
        
        return vmanager


class WorkflowVariableManagerInitialiser(VariableManagerInitialiser):
    
    def initialise(self) -> VariableManager:
        vmanager = VariableManager()
        
        for tinput in self.tool.tool_inputs():
            if not task_inputs.exists(self.tool.id(), tinput):
                self.initialise_ignored_input(tinput, vmanager)
                continue
            
            task_input = task_inputs.get(self.tool.id(), tinput)
            if task_input.ti_type == TaskInputType.IGNORED:
                self.initialise_ignored_input(tinput, vmanager)
            else:
                vtype_str = self.types_map[task_input.ti_type]
                vmanager.update(task_input.tinput_id, vtype_str=vtype_str, value=task_input.value)
        
        return vmanager
    
    def initialise_ignored_input(self, tinput: TInput, vmanager: VariableManager) -> None:
        if tinput.default is not None:
            vtype_str = 'static'
            vmanager.update(tinput.id(), vtype_str=vtype_str, value=tinput.default)
        else:
            vtype_str = 'ignored'
            vmanager.update(tinput.id(), vtype_str=vtype_str, value=None)


class VariableType(Enum): 
    TASK_INPUT  = auto()
    PARAM       = auto()
    STATIC      = auto()
    IGNORED     = auto()
    CHANNEL     = auto()
    LOCAL       = auto()

@dataclass
class Variable:
    vtype: VariableType
    value: Any

@dataclass
class VariableHistory:
    tinput_id: str
    items: list[Variable] = field(default_factory=list)

    @property
    def original(self) -> Variable:
        return self.items[0]
    
    @property
    def previous(self) -> Variable:
        return self.items[-2]
    
    @property
    def current(self) -> Variable:
        return self.items[-1]
    


class VariableManager:
    """
    within a process, for a given TInput, get the variable which currently represents that TInput.
    eg 

    input:
    1 path indexed_bam_flat
    2 path cath_crossmapped, stageAs: 'cath_crossmapped'

    script:
    3 def in_bams = get_primary_files(indexed_bam_flat, 2)
    4 def in_bams_joined = in_bams.join(' ')
    5 def cath_crossmapped = cath_crossmapped ? "-cx ${cath_crossmapped}" : ""
    '''
    6 echo \
    7 ${in_bams_joined} \
    8 ${cath_crossmapped} \
    9 --promoter=${params.promoter_bp} \
    '''

    TInput('inBams')
    line 1 -> indexed_bam_flat
    line 3 -> in_bams
    line 4 -> in_bams_joined
    line 6 -> in_bams_joined

    TInput('cathCrossmapped')
    line 2 -> cath_crossmapped
    line 8 -> cath_crossmapped

    TInput('promoterBp')
    line 9 -> params.promoter_bp

    data_structure description:
    {
        tinput_id: [
            first,  (original)
            second, 
            third,  (current)
            ..
        ]
    }

    """
    def __init__(self) -> None:
        self.data_structure: dict[str, VariableHistory] = {}
        self.type_map = {
            'task_input': VariableType.TASK_INPUT,
            'param': VariableType.PARAM,
            'static': VariableType.STATIC,
            'ignored': VariableType.IGNORED,
            'channel': VariableType.CHANNEL,
            'local': VariableType.LOCAL,
        }

    def get(self, tinput_id: str) -> VariableHistory:
        if tinput_id not in self.data_structure:
            raise RuntimeError
        return self.data_structure[tinput_id]
        
    def update(self, tinput_id: str, vtype_str: str, value: Optional[str | list[str]]) -> None:
        # create new variable history to track tinput_id if not exists
        if tinput_id not in self.data_structure:
            new_history = VariableHistory(tinput_id)
            self.data_structure[tinput_id] = new_history
        
        # gen new variable & update history
        history = self.data_structure[tinput_id]
        vtype = self.type_map[vtype_str]
        new_variable = Variable(vtype, value)
        history.items.append(new_variable)

    def to_string(self) -> str:
        out: str = ''
        for tinput_id, history in self.data_structure.items():
            block = ''
            block += f'{tinput_id}:\n'
            for var in history.items:
                block += f'[{var.vtype}] {var.value}\n'
            block += '\n'
            out += block
        return out
    
