

from janis_core import WorkflowBuilder

from .gather import gather_tools
from .tools import prune_tools_and_sources
from .workflows import prune_main_workflow_inputs


"""
exists to enable the settings.translate.MODE feature. 
prunes the main workflow inputs, so that only those which are used are kept.
prunes tool inputs which aren't used. 
prunes step sources which aren't used.

"""
 
def prune_workflow(wf: WorkflowBuilder) -> None:
    tools = gather_tools(wf)
    prune_tools_and_sources(wf, tools)
    prune_main_workflow_inputs(wf)


