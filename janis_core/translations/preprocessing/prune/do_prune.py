

from janis_core import Workflow, Tool, ToolInput, CommandToolBuilder

from .history import TaskInputCollector
from ...common import trace

def prune_unused_tool_inputs(wf: Workflow, tools: dict[str, Tool]) -> dict[str, Tool]:
    pruned_tools: dict[str, Tool] = {}
    for tool_id, tool in tools.items():
        pruned_tools[tool_id] = prune(wf, tool)
    return pruned_tools

def prune(wf: Workflow, tool: Tool) -> Tool:
    if not isinstance(tool, CommandToolBuilder):
        return tool
    
    # get the tinputs which are not needed based on step inputs
    collector = TaskInputCollector(tool)
    collector.collect(wf)
    step_input_ids = get_used_step_tinputs(collector)

    # get the tinputs which are referenced in tool outputs
    output_referenced_ids = get_output_referenced_tinputs(tool)

    # remove the tinputs from the tool using this data
    tool = prune_tool(tool, step_input_ids, output_referenced_ids)

    return tool

def get_used_step_tinputs(collector: TaskInputCollector) -> set[str]:
    tinput_to_keep: set[str] = set()
    
    for tinput_id, history in collector.histories.items():
        # RULE 1: mandatory tool inputs are always kept
        if not history.is_optional:  
            tinput_to_keep.add(tinput_id)
        # RULE 2: anything which is supplied a non-null value must be kept
        elif len(history.non_null_unique_values) > 0:
            tinput_to_keep.add(tinput_id)
        else:
            continue
    
    return tinput_to_keep

def get_output_referenced_tinputs(tool: CommandToolBuilder) -> set[str]:
    referenced_tinputs: set[str] = set()
    for tout in tool._outputs:
        referenced_tinputs = referenced_tinputs | trace.trace_referenced_variables(tout, tool)
    return referenced_tinputs

def prune_tool(tool: CommandToolBuilder, step_input_ids: set[str], output_referenced_ids: set[str]) -> Tool:
    new_inputs: list[ToolInput] = []
    for tinput in tool._inputs:     # type: ignore
        if tinput.id() in step_input_ids or tinput.id() in output_referenced_ids:
            new_inputs.append(tinput)
    tool._inputs = new_inputs       # type: ignore
    return tool

