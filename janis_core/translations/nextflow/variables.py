
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Optional
from enum import Enum, auto
from dataclasses import dataclass, field

from janis_core import Tool, Workflow, CommandTool, PythonTool, TInput

from . import task_inputs
from .task_inputs import TaskInputType




def init_variable_manager_for_task(tool: Tool) -> VariableManager:
    """
    initialise a variable manager for a tool, covering all tool inputs. 
    differs from task_inputs because task_inputs are based on the task 'class', whereas the variable manager is specific to an 'instance' of the task 
    """
    if isinstance(tool, Workflow):
        initialiser = WorkflowVariableManagerInitialiser(tool)
    else:
        initialiser = ProcessVariableManagerInitialiser(tool)
    return initialiser.initialise()



    # # process & workflow
    # all_tinputs = set(tool.tool_inputs())
    # all_task_inputs = set([x.tinput_id for x in task_inputs.getall(tool.id())])
    
    # for tinput in all_tinputs:
    #     # if the tinput is a task_input (ie constant per instance), add it to the variable manager
    #     if tinput.id() in all_task_inputs:
    #         task_input = task_inputs.get(tool.id(), tinput)
    #         vtype_str = types_map[task_input.ti_type]
    #         vmanager.update(task_input.tinput_id, vtype_str=vtype_str, value=task_input.value)
        
    #     # process instance
    #     elif isinstance(tool, CommandTool | PythonTool):
    #         # i dont think this path is ever taken
    #         vmanager.update(tinput.id(), vtype_str='ignored', value=None)

    #     # workflow instance
    #     else:
    #         if tinput.default is not None:
    #             vmanager.update(tinput.id(), vtype_str='static', value=tinput.default)
    #         else:
    #             pass
        
    # return vmanager


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
            
            else:
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


    

    # param_tinput_ids = get_true_workflow_inputs(wf)
    # static_tinput_ids = get_static_workflow_inputs(all_tinput_ids, param_tinput_ids, wf)
    # ignored_tinput_ids = get_ignored_workflow_inputs(all_tinput_ids, param_tinput_ids, static_tinput_ids)
    
    # # param inputs
    # for tinput_id in param_tinput_ids:
    #     ti_type = 'param'
    #     tinput = [x for x in wf.tool_inputs() if x.id() == tinput_id][0]
    #     param = params.register(tinput, task_id=wf.id())
    #     value = f'params.{param.name}'
    #     task_inputs.update(wf.id(), ti_type, tinput_id, value)
    #     print()
    
    # # static inputs
    # for tinput_id in static_tinput_ids:
    #     ti_type = 'static'
    #     tinput = [x for x in wf.tool_inputs() if x.id() == tinput_id][0]
    #     value = tinput.default
    #     task_inputs.update(wf.id(), ti_type, tinput_id, value)
    #     print()
    
    # # ignored inputs
    # for tinput_id in ignored_tinput_ids:
    #     ti_type = 'ignored'
    #     value = None
    #     task_inputs.update(wf.id(), ti_type, tinput_id, value)
    #     print()




# def get_static_workflow_inputs(all_tinput_ids: set[str], param_tinput_ids: set[str], wf: Workflow) -> set[str]:
#     out: set[str] = set()
#     for tinput_id in all_tinput_ids:
#         if tinput_id not in param_tinput_ids:
#             tinput = [x for x in wf.tool_inputs() if x.id() == tinput_id][0]
#             if tinput.default is not None:
#                 out.add(tinput_id)
#     return out

# def get_ignored_workflow_inputs(all_tinput_ids: set[str], param_tinput_ids: set[str], static_tinput_ids: set[str]) -> set[str]:
#     surviving_ids = deepcopy(all_tinput_ids)
#     surviving_ids -= param_tinput_ids
#     surviving_ids -= static_tinput_ids
#     return surviving_ids
    


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
    value: Optional[str | list[str]]

@dataclass
class VariableHistory:
    tinput_id: str
    items: list[Variable] = field(default_factory=list)

    @property
    def original(self) -> Variable:
        return self.items[0]
    
    @property
    def current(self) -> Variable:
        return self.items[-1]
    
    @property
    def all(self) -> list[Variable]:
        return self.items


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
    
