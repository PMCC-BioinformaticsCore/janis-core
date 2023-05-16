

from janis_core import Workflow, Tool, ToolInput, CommandToolBuilder

from .history import TaskInputCollector


def prune_unused_tool_inputs(wf: Workflow, tools: dict[str, Tool]) -> dict[str, Tool]:
    pruned_tools: dict[str, Tool] = {}
    for tool_id, tool in tools.items():
        pruned_tools[tool_id] = prune(wf, tool)
    return pruned_tools

def prune(wf: Workflow, tool: Tool) -> Tool:
    collector = TaskInputCollector(tool)
    collector.collect(wf)

    tinputs_to_remove = get_tinputs_to_remove(collector)
    if isinstance(tool, Workflow):
        pass
        # tool = prune_workflow(tool, tinputs_to_remove)
    elif isinstance(tool, CommandToolBuilder):
        tool = prune_tool(tool, tinputs_to_remove)
    return tool

def get_tinputs_to_remove(collector: TaskInputCollector) -> set[str]:
    tinputs_to_remove: set[str] = set()
    
    for tinput_id, history in collector.histories.items():
        # RULE 1: mandatory tool inputs are always kept
        if not history.is_optional:  
            continue
        # RULE 2: anything which is supplied a non-null value must be kept
        elif len(history.non_null_unique_values) > 0:
            continue
        else:
            tinputs_to_remove.add(tinput_id)
    
    return tinputs_to_remove

# def prune_workflow(wf: Workflow, tinputs_to_remove: set[str]) -> Workflow:
#     for tinput_id in tinputs_to_remove:
#         del wf.input_nodes[tinput_id]
#         del wf.nodes[tinput_id]
#     return wf

def prune_tool(tool: CommandToolBuilder, tinputs_to_remove: set[str]) -> Tool:
    new_inputs: list[ToolInput] = []
    for tinput in tool._inputs:     # type: ignore
        if tinput.id() not in tinputs_to_remove:
            new_inputs.append(tinput)
    tool._inputs = new_inputs       # type: ignore
    return tool

