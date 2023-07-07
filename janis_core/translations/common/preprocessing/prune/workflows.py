

from janis_core import WorkflowBuilder
from janis_core.operators import InputNodeSelector
from janis_core.translations.common import trace


def prune_main_workflow_inputs(wf: WorkflowBuilder) -> None:
    valid_workflow_inputs: set[str] = set()
    valid_workflow_inputs = valid_workflow_inputs | get_mandatory_input_ids(wf)
    valid_workflow_inputs = valid_workflow_inputs | get_referenced_input_ids(wf)
    invalid_workflow_inputs = set(wf.input_nodes.keys()) - valid_workflow_inputs
    for tinput_id in invalid_workflow_inputs:
        del wf.input_nodes[tinput_id]

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
            entities = trace.trace_entities(src)
            for entity in entities:
                if isinstance(entity, InputNodeSelector):
                    collected_ids.add(entity.input_node.id())

    return collected_ids