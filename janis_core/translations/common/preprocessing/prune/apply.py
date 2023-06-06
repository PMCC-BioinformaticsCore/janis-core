

from janis_core import Workflow
from janis_core import Tool


def apply_pruned_tools(wf: Workflow, pruned_tools: dict[str, Tool]) -> Workflow:
    for step in wf.step_nodes.values():
        if isinstance(step.tool, Workflow):
            step.tool = apply_pruned_tools(step.tool, pruned_tools)
        elif step.tool.id() in pruned_tools:
            step.tool = pruned_tools[step.tool.id()]
        else:
            continue
    return wf
