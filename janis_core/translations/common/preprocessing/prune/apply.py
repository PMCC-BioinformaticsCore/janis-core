

from janis_core import WorkflowBase
from janis_core import Tool


def apply_pruned_tools(wf: WorkflowBase, pruned_tools: dict[str, Tool]) -> WorkflowBase:
    for step in wf.step_nodes.values():
        if isinstance(step.tool, WorkflowBase):
            step.tool = apply_pruned_tools(step.tool, pruned_tools)
        elif step.tool.id() in pruned_tools:
            step.tool = pruned_tools[step.tool.id()]
        else:
            continue
    return wf
