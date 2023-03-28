from typing import Any
from copy import deepcopy
from dataclasses import dataclass
from abc import ABC, abstractmethod

from janis_core import Workflow

from ..scope import Scope
from ..plumbing import gen_task_call

from .. import naming
from .. import channels
from .. import unwrap

from ..model.process import NFProcess
from ..model.files import NFImportsBlock
from ..model.files import NFFunctionsBlock
from ..model.workflow import NFWorkflow
from ..model.workflow import NFMainWorkflow
from ..model.workflow import NFSubWorkflow
from ..model.workflow import NFWorkflowTake
from ..model.workflow import NFWorkflowEmit



def gen_workflow(name: str, scope: Scope, wf: Workflow, item_register: Any) -> NFWorkflow:
    """Generate a Nextflow Workflow object"""
    is_subworkflow = True if len(scope.labels) > 1 else False
    if is_subworkflow:
        generator = MainWFGenerator(scope, name, wf, item_register)
    else:
        generator = SubWFGenerator(scope, name, wf, item_register)
    return generator.generate()


@dataclass 
class WFGenerator(ABC):
    scope: Scope
    name: str
    wf: Workflow
    item_register: Any

    @property
    def nf_name(self) -> str:
        return naming.constructs.gen_varname_workflow(self.name)
    
    @property
    def main_block(self) -> list[str]:
        # MAIN (workflow step calls, channel operations)
        main: list[str] = []
        for step in self.wf.step_nodes.values():
            current_scope = deepcopy(self.scope)
            current_scope.update(step)
            nf_items = self.item_register.get(current_scope)
            
            for nf_item in nf_items:
                if isinstance(nf_item, channels.ChannelOperation):
                    main.append(nf_item.get_string())
                    continue
                elif isinstance(nf_item, NFProcess) or isinstance(nf_item, NFWorkflow):
                    entity_name = nf_item.name
                elif isinstance(nf_item, NFImportsBlock):
                    continue
                elif isinstance(nf_item, NFFunctionsBlock):
                    continue
                # future items
                else:
                    raise NotImplementedError
            
                task_call = gen_task_call(step, current_scope, entity_name)
                main.append(task_call)
        return main
    
    @abstractmethod
    def generate(self) -> NFWorkflow:
        ...


@dataclass 
class MainWFGenerator(WFGenerator):

    def generate(self) -> NFMainWorkflow:
        return NFMainWorkflow(self.nf_name, self.main_block)


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
        return NFSubWorkflow(self.nf_name, self.main_block, self.take_items, self.emit_items)
   
