

from janis_core import WorkflowBase
from janis_core import Tool


def gather_tools(wf: WorkflowBase) -> dict[str, Tool]:
    gathered_tools: dict[str, Tool] = {}
    return do_gather_tools(wf, gathered_tools)

def do_gather_tools(wf: WorkflowBase, gathered_tools: dict[str, Tool]) -> dict[str, Tool]:
    for step in wf.step_nodes.values():
        if isinstance(step.tool, WorkflowBase):
            # if step.tool.id() not in gathered_tools:
            #     gathered_tools[step.tool.id()] = step.tool
            gathered_tools = do_gather_tools(step.tool, gathered_tools)
        elif step.tool.id() not in gathered_tools:
            gathered_tools[step.tool.id()] = step.tool
        else:
            continue
    return gathered_tools

