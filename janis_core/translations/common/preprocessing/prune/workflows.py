

from janis_core import WorkflowBuilder
from janis_core.operators import InputNodeSelector

def prune_main_workflow_inputs(wf: WorkflowBuilder) -> WorkflowBuilder:
    valid_workflow_inputs: set[str] = set()

    for step in wf.step_nodes.values():
        for src in step.sources.values():
            selector = src.source_map[0].source
            if isinstance(selector, InputNodeSelector):
                valid_workflow_inputs.add(selector.input_node.id())
    
    wf.input_nodes = {k: v for k, v in wf.input_nodes.items() if k in valid_workflow_inputs}
    return wf

