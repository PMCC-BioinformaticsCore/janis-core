

from janis_core import WorkflowBase

from .gather import gather_tools
from .do_prune import prune_unused_tool_inputs
from .apply import apply_pruned_tools

 
def prune_unused_inputs(wf: WorkflowBase) -> WorkflowBase:
    tools = gather_tools(wf)
    pruned_tools = prune_unused_tool_inputs(wf, tools)
    wf = apply_pruned_tools(wf, pruned_tools)
    return wf
