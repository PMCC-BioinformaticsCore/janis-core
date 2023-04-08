
from copy import deepcopy
from dataclasses import dataclass
from abc import ABC, abstractmethod

from janis_core import Workflow
from janis_core import settings

from .. import naming
from .. import channels
from .. import unwrap

from ..plumbing import gen_task_call
from ..casefmt import to_case
from ..variables import init_variable_manager_for_task

from ..model.process import NFProcess
from ..model.files import NFImportsBlock
from ..model.files import NFFunctionsBlock
from ..model.workflow import NFWorkflow
from ..model.workflow import NFMainWorkflow
from ..model.workflow import NFSubWorkflow
from ..model.workflow import NFWorkflowTake
from ..model.workflow import NFWorkflowEmit


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
        assert(tool_id) not in self.workflow_dict
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
    
    @property
    def name(self) -> str:
        return naming.constructs.gen_varname_workflow(self.wf.id())

    def __post_init__(self) -> None:
        self.vmanager = init_variable_manager_for_task(self.wf)
        print()

    @property
    def main_block(self) -> list[str]:
        # MAIN (workflow step calls, channel operations)
        main: list[str] = []
        for step in self.wf.step_nodes.values():
            # workflow_scope is where we are now
            workflow_scope = self.scope
            
            # task_scope is the scope of the subtask 
            task_scope = deepcopy(self.scope)
            task_scope.update(step)
            nf_items = self.item_register.get(task_scope)
            
            # ordered items which appear in the workflow block (between workflow {...})
            # only interested in channel operations and tasks to call
            for nf_item in nf_items:

                # channel operation
                if isinstance(nf_item, channels.ChannelOperation):
                    main.append(nf_item.get_string())
                
                # subtask to call
                elif isinstance(nf_item, NFProcess) or isinstance(nf_item, NFWorkflow):
                    entity_name = to_case(step.id(), settings.translate.nextflow.NF_PROCESS_CASE)
                    task_call = gen_task_call(step, workflow_scope, entity_name)
                    main.append(task_call)
                
                elif isinstance(nf_item, NFImportsBlock):
                    continue
                
                elif isinstance(nf_item, NFFunctionsBlock):
                    continue
                
                # future items
                else:
                    raise NotImplementedError

        return main
    
    @abstractmethod
    def generate(self) -> NFWorkflow:
        ...


@dataclass 
class MainWFGenerator(WFGenerator):

    def generate(self) -> NFMainWorkflow:
        return NFMainWorkflow(self.name, self.alias, self.main_block)


@dataclass 
class SubWFGenerator(WFGenerator):

    @property
    def take_items(self) -> list[NFWorkflowTake]:
        relevant_channels = channels.getall(self.scope)
        return [NFWorkflowTake(ch.name) for ch in relevant_channels]
    
    @property
    def emit_items(self) -> list[NFWorkflowEmit]:
        emit: list[NFWorkflowEmit] = []
        for out in self.wf.output_nodes.values():
            outname = out.id()
            expression = unwrap.unwrap_expression(
                val=out.source, 
                context='workflow',
                scope=self.scope, 
                in_shell_script=True
            )
            emit.append(NFWorkflowEmit(outname, expression))
        return emit

    def generate(self) -> NFSubWorkflow:
        return NFSubWorkflow(self.name, self.alias, self.main_block, self.take_items, self.emit_items)
   
