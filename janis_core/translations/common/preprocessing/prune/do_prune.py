

from janis_core import WorkflowBuilder, Tool, ToolInput, CommandToolBuilder,ToolArgument

from .history import TaskInputCollector
from janis_core.translations.common import trace


def prune_unused_tool_inputs(wf: WorkflowBuilder, tools: dict[str, Tool]) -> dict[str, Tool]:
    pruned_tools: dict[str, Tool] = {}
    for tool_id, tool in tools.items():
        pruned_tools[tool_id] = prune(wf, tool)
    return pruned_tools

def prune(wf: WorkflowBuilder, tool: Tool) -> Tool:
    if not isinstance(tool, CommandToolBuilder):
        return tool
    
    valid_tinput_ids: set[str] = set()

    # get the tinputs which are referenced in step inputs
    valid_tinput_ids = valid_tinput_ids | get_used_step_tinputs(wf, tool)
    
    # get the tinputs which are referenced in tool outputs
    valid_tinput_ids = valid_tinput_ids | get_output_referenced_tinputs(wf, tool)
    
    # get the tinputs which reference previously validified tinputs
    valid_tinput_ids = valid_tinput_ids | get_tinput_reference_tinputs(tool, valid_tinput_ids)

    # remove the tinputs from the tool using this data
    tool = prune_tool(tool, valid_tinput_ids)

    return tool

def get_used_step_tinputs(wf: WorkflowBuilder, tool: CommandToolBuilder) -> set[str]:
    # get the tinputs which are not needed based on step inputs
    collector = TaskInputCollector(tool)
    collector.collect(wf)
    tinput_to_keep: set[str] = set()
    
    for tinput_id, history in collector.histories.items():
        # RULE 1: mandatory tool inputs are always kept
        if not history.is_optional:  
            tinput_to_keep.add(tinput_id)
        # RULE 2: step connections are always kept
        elif history.connections:
            tinput_to_keep.add(tinput_id)
        # RULE 3: anything which is supplied a non-null value must be kept
        elif len(history.genuine_input_sources) > 0:
            tinput_to_keep.add(tinput_id)
        else:
            continue
    
    return tinput_to_keep

def get_tinput_reference_tinputs(tool: CommandToolBuilder, valid_tinput_ids: set[str]) -> set[str]:
    """
    Gets the tinputs which have a reference to step input tinputs
    eg: 
        step_input_ids = {'myfile'}
        tool inputs = [
        ToolInput("myfile", File()),
        ToolInput("myfilename", Filename(prefix=InputSelector("myfile"), extension=""))
        ]
    we should keep "myfilename" because it references "myfile"
    """
    extra_tinput_ids: set[str] = set()
    for tinput in tool._inputs:
        # early exit for previously validated tinputs
        if tinput.id() in valid_tinput_ids:
            continue
        
        # reference tracing for unvalidated tinputs - do they link to a validated tinput?
        refs = trace.trace_referenced_variables(tinput, tool)
        for ref in refs:
            if ref in valid_tinput_ids:
                extra_tinput_ids.add(tinput.id())
                break
    return extra_tinput_ids

def get_output_referenced_tinputs(wf: WorkflowBuilder, tool: CommandToolBuilder) -> set[str]:
    extra_tinput_ids: set[str] = set()
    for tout in tool._outputs:
        extra_tinput_ids = extra_tinput_ids | trace.trace_referenced_variables(tout, tool)
    return extra_tinput_ids

def prune_tool(tool: CommandToolBuilder, valid_tinput_ids: set[str]) -> Tool:
    new_inputs = prune_tool_inputs(tool, valid_tinput_ids)       
    new_arguments = prune_tool_arguments(tool, valid_tinput_ids) 
    tool._inputs = new_inputs        # type: ignore
    tool._arguments = new_arguments  # type: ignore
    return tool
    
def prune_tool_inputs(tool: CommandToolBuilder, valid_tinput_ids: set[str]) -> list[ToolInput]:
    new_inputs: list[ToolInput] = []
    if tool._inputs is None:        # type: ignore
        return new_inputs
    for tinput in tool._inputs:     # type: ignore
        if tinput.id() in valid_tinput_ids:
            new_inputs.append(tinput)
    return new_inputs

def prune_tool_arguments(tool: CommandToolBuilder, valid_tinput_ids: set[str]) -> list[ToolArgument]:
    new_args: list[ToolArgument] = []
    if tool._arguments is None:     # type: ignore
        return new_args
    for targ in tool._arguments:     # type: ignore
        refs = trace.trace_referenced_variables(targ, tool)
        if all([ref in valid_tinput_ids for ref in refs]):
            new_args.append(targ)
    return new_args
    