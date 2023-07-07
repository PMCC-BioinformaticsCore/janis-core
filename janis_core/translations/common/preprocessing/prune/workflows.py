

from dataclasses import dataclass, field

from janis_core import WorkflowBuilder
from janis_core.types import DataType
from janis_core.workflow.workflow import InputNode
from janis_core.operators import InputNodeSelector
from janis_core.translations.common import trace


def prune_main_workflow_inputs(wf: WorkflowBuilder) -> None:
    valid_workflow_inputs: set[str] = set()
    valid_workflow_inputs = valid_workflow_inputs | get_mandatory_input_ids(wf)
    valid_workflow_inputs = valid_workflow_inputs | get_referenced_input_ids(wf)
    wf.input_nodes = {k: v for k, v in wf.input_nodes.items() if k in valid_workflow_inputs}

@dataclass 
class InputNodeReference:
    step_id: str 
    tinput_id: str
    tinput_dtype: DataType

@dataclass 
class InputNodeReferenceCatalogue:
    node: InputNode
    references: list[InputNodeReference] = field(default_factory=list)
    
    def add(self, step_id: str, tinput_id: str, tinput_dtype: DataType):
        self.references.append(InputNodeReference(step_id, tinput_id, tinput_dtype))

def get_mandatory_input_ids(wf: WorkflowBuilder) -> set[str]:
    collected_ids: set[str] = set()
    
    for tinput_id, tinput in wf.input_nodes.items():
        if tinput.datatype.optional == False:
            collected_ids.add(tinput_id)
    
    return collected_ids

def get_referenced_input_ids(wf: WorkflowBuilder) -> set[str]:
    collected_ids: set[str] = set()
    
    for step in wf.step_nodes.values():
        for src in step.sources.values():
            selector = src.source_map[0].source
            if isinstance(selector, InputNodeSelector):
                collected_ids.add(selector.input_node.id())
    
    return collected_ids
