
from dataclasses import dataclass
from abc import ABC, abstractmethod

from janis_core import settings
from janis_core import TInput
from janis_core.workflow.workflow import StepNode, Workflow
from janis_core.types import File

from ... import naming
from ... import unwrap

from ...casefmt import to_case
from ...variables import init_variable_manager_for_task
from ...variables import VariableType

from ...model.process import NFProcess
from ...model.workflow import NFWorkflow
from ...model.workflow import NFMainWorkflow
from ...model.workflow import NFSubWorkflow
from ...model.workflow import NFWorkflowTake
from ...model.workflow import NFWorkflowEmit

from janis_core import translation_utils as utils

from .call import gen_task_call


def generate_workflows(wf: Workflow, process_dict: dict[str, NFProcess]) -> dict[str, NFWorkflow]:
    """
    for each workflow in workflow, generate a nextflow workflow.
    includes subworkflows, & main workflow.
    """
    manager = WorkflowGenerationManager(wf, process_dict)
    return manager.generate()


class WorkflowGenerationManager:
    def __init__(self, wf: Workflow, process_dict: dict[str, NFProcess]) -> None:
        self.main_wf = wf
        self.process_dict = process_dict
        self.workflow_dict: dict[str, NFWorkflow] = {}
    
    def generate(self) -> dict[str, NFWorkflow]:
        self.do_generate(self.main_wf, is_subworkflow=False)
        return self.workflow_dict

    def do_generate(self, wf: Workflow, is_subworkflow: bool=True) -> None:
        """
        handles depth-first recursive workflow parsing. 
        subworkflows are parsed before calling workflows, so that the subworkflow can be called. 
        """
        is_main = True if not is_subworkflow else False
        is_subworkflow = True

        # depth first 
        for step in wf.step_nodes.values():
            if isinstance(step.tool, Workflow):
                # recursively do for sub-sub-subworkflows 
                self.do_generate(step.tool, is_subworkflow=True)
                # do for sub-subworkflow
                tool_id = step.tool.id()
                if tool_id not in self.workflow_dict: 
                    workflow = self.generate_workflow(step.tool, is_subworkflow=True)
                    self.workflow_dict[tool_id] = workflow
        
        # do for this subworkflow
        tool_id = wf.id()
        if tool_id not in self.workflow_dict:
            workflow = self.generate_workflow(wf, is_subworkflow = not is_main)  # im sorry for this, I hate it too. 
            self.workflow_dict[tool_id] = workflow
    
    def generate_workflow(self, wf: Workflow, is_subworkflow: bool=False) -> NFWorkflow:
        """Generate a Nextflow Workflow object"""
        if is_subworkflow:
            generator = SubWFGenerator(wf, self.workflow_dict, self.process_dict)
        else:
            generator = MainWFGenerator(wf, self.workflow_dict, self.process_dict)
        return generator.generate()


@dataclass 
class WFGenerator(ABC):
    wf: Workflow
    workflow_dict: dict[str, NFWorkflow]
    process_dict: dict[str, NFProcess]

    def __post_init__(self) -> None:
        self.vmanager = init_variable_manager_for_task(self.wf)
        print()
    
    @property
    def name(self) -> str:
        return naming.constructs.gen_varname_workflow(self.wf.id())
    
    @property
    @abstractmethod
    def take_block(self) -> list[NFWorkflowTake]:
        ...
    
    @property
    @abstractmethod
    def emit_block(self) -> list[NFWorkflowEmit]:
        ...
    
    @property
    def main_block(self) -> list[str]:
        # MAIN (workflow step calls, channel operations)
        main: list[str] = []
        for step in self.wf.step_nodes.values():
            main += self.gen_operations(step)
            main += self.gen_call(step)
        return main

    def gen_operations(self, step: StepNode) -> list[str]:
        # in case any channel operations need to occur before task 
        # call to ensure correct plumbing. 
        # in these cases, the variable name used to reference the tinput 
        # will be updated in self.vmanager
        # nothing for now though. 
        out: list[str] = []
        return out
    
    def gen_call(self, step: StepNode) -> list[str]:
        call: list[str] = []
    
        task_id = step.tool.id()
        alias = to_case(step.id(), settings.translate.nextflow.NF_PROCESS_CASE)

        if task_id in self.process_dict:
            task = self.process_dict[task_id]
        elif task_id in self.workflow_dict:
            task = self.workflow_dict[task_id]
        else:
            raise RuntimeError('task has not been generated')
        
        task_call = gen_task_call(alias, task, self.vmanager, step)
        call += task_call
        call += ['']
    
        return call

    @abstractmethod
    def generate(self) -> NFWorkflow:
        ...


@dataclass 
class MainWFGenerator(WFGenerator):

    @property
    def take_block(self) -> list[NFWorkflowTake]:
        return []
    
    @property
    def emit_block(self) -> list[NFWorkflowEmit]:
        return []
    
    # @property
    # def param_inputs(self) -> set[str]:
    #     out: set[str] = set()
    #     for tinput in self.wf.tool_inputs():
    #         var = self.vmanager.get(tinput.id()).current
    #         if var.vtype == VariableType.PARAM:
    #             out.add(tinput.id())
    #     return out
    
    @property
    def param_inputs(self) -> list[TInput]:
        out: list[TInput] = []
        for tinput in self.wf.tool_inputs():
            var = self.vmanager.get(tinput.id()).current
            if var.vtype == VariableType.PARAM:
                out.append(tinput)
        return out
    
    def update_variables(self) -> None:
        for tinput in self.param_inputs:
            if utils.is_file_type(tinput.intype):
                if tinput.intype.optional:
                    f_name = naming.constructs.gen_varname_file(tinput.id(), dtype=tinput.intype)
                    self.vmanager.update(tinput.id(), 'local', f_name)
                else:
                    ch_name = naming.constructs.gen_varname_channel(tinput.id(), dtype=tinput.intype)
                    self.vmanager.update(tinput.id(), 'channel', ch_name)

    def generate(self) -> NFMainWorkflow:
        self.update_variables()
        return NFMainWorkflow(self.name, self.main_block, self.take_block, self.emit_block)


@dataclass 
class SubWFGenerator(WFGenerator):

    @property
    def take_ids(self) -> set[str]:
        out: set[str] = set()
        for tinput in self.wf.tool_inputs():
            cvar = self.vmanager.get(tinput.id()).current
            if cvar.vtype == VariableType.TASK_INPUT:
                out.add(tinput.id())
        return out

    @property
    def take_block(self) -> list[NFWorkflowTake]:
        take: list[NFWorkflowTake] = []
        for tinput_id in self.take_ids:
            name = naming.constructs.gen_varname_channel(tinput_id)
            dtype = [x for x in self.wf.tool_inputs() if x.id() == tinput_id][0].intype  # type: ignore
            take_item = NFWorkflowTake(name, tinput_id, dtype)
            take.append(take_item)
        return take
    
    @property
    def emit_block(self) -> list[NFWorkflowEmit]:
        emit: list[NFWorkflowEmit] = []
        for out in self.wf.output_nodes.values():
            outname = out.id()
            expression = unwrap.unwrap_expression(
                val=out.source, 
                context='workflow',
                variable_manager=self.vmanager,
                in_shell_script=True
            )
            emit_item = NFWorkflowEmit(outname, expression)
            emit.append(emit_item)
        return emit

    def generate(self) -> NFSubWorkflow:
        return NFSubWorkflow(self.name, self.main_block, self.take_block, self.emit_block)
   
