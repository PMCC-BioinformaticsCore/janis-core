


from janis_core import WorkflowBuilder, CommandToolBuilder, ToolInput, ToolArgument
from janis_core.operators.selectors import InputNodeSelector
from janis_core.translations.common import trace

from .history import TaskInputCollector


def prune_tools_and_sources(wf: WorkflowBuilder, tools: dict[str, CommandToolBuilder]) -> WorkflowBuilder:
    for tool_id, tool in tools.items():
        wf = prune(wf, tool)
    return wf

def prune(wf: WorkflowBuilder, tool: CommandToolBuilder) -> WorkflowBuilder:
    # for each tool input, get all sources from each wf step which feed this input 
    collector = TaskInputCollector(tool)
    collector.collect(wf)

    # identifying which tinputs need to be fed via step.sources
    valid_tinput_ids: set[str] = set() 
    # identify tinputs which need to be kept based on step.sources
    valid_tinput_ids = valid_tinput_ids | get_step_referenced_tinputs(collector)
    # identify tinputs which are referenced in tool outputs
    valid_tinput_ids = valid_tinput_ids | get_output_referenced_tinputs(wf, tool)

    # for tinputs which only have a single static value as source, transfer the static value to tool default
    tool = migrate_statics_to_defaults(tool, collector, valid_tinput_ids)
    # remove step.sources which are not needed
    wf = prune_sources(wf, tool, valid_tinput_ids)

    # add tinputs which reference a previously validated tinput
    valid_tinput_ids = valid_tinput_ids | get_tinput_reference_tinputs(tool, valid_tinput_ids)
    # add tinputs which have a non-null default
    valid_tinput_ids = valid_tinput_ids | get_default_tinputs(tool, valid_tinput_ids)
    
    # remove the ToolInputs / ToolArguments which are not needed
    tool = prune_tool(tool, valid_tinput_ids)
    # apply the pruned tool to the workflow
    wf = apply_pruned_tool(wf, tool)

    return wf

def prune_sources(wf: WorkflowBuilder, tool: CommandToolBuilder, valid_tinput_ids: set[str]) -> WorkflowBuilder:
    for step in wf.step_nodes.values():
        if isinstance(step.tool, WorkflowBuilder):
            step.tool = prune_sources(step.tool, tool, valid_tinput_ids)
        elif step.tool.id() == tool.id():
            step.sources = {k: v for k, v in step.sources.items() if k in valid_tinput_ids}
        else:
            continue
    return wf

def get_step_referenced_tinputs(collector: TaskInputCollector) -> set[str]:
    # get the tinputs which are needed based on step inputs
    tinputs_to_keep: set[str] = set()
    
    for tinput_id, history in collector.histories.items():
        # RULE 1: mandatory tool inputs are always kept
        if not history.is_optional:  
            tinputs_to_keep.add(tinput_id)
        # RULE 2: if has step connections, keep
        elif history.connections:
            tinputs_to_keep.add(tinput_id)
        # RULE 3: if has 1+ sources which are mandatory workflow inputs, keep
        elif len(history.mandatory_input_sources) >= 1:
            tinputs_to_keep.add(tinput_id)
        # RULE 4: if has multiple sources, keep
        elif len(history.sources) >= 2:
            tinputs_to_keep.add(tinput_id)
        # RULE 5: (edge case) if tool used in 2+ steps, but tinput has only 1 source, keep.
        #         this is needed because if we are driving the tinput's value from a source in 
        #         one step, then it's value is driven by its default in another. they could be different. 
        elif len(history.sources) == 1 and collector.step_count >= 2:
            tinputs_to_keep.add(tinput_id)
        else:
            continue
    
    return tinputs_to_keep

def get_output_referenced_tinputs(wf: WorkflowBuilder, tool: CommandToolBuilder) -> set[str]:
    # TODO check if the output is consumed as a source somewhere else in the workflow. 
    # if not we can remove it. 
    extra_tinput_ids: set[str] = set()
    for tout in tool._outputs:
        extra_tinput_ids = extra_tinput_ids | trace.trace_referenced_variables(tout, tool)
    return extra_tinput_ids

def migrate_statics_to_defaults(
    tool: CommandToolBuilder, 
    collector: TaskInputCollector, 
    valid_tinput_ids: set[str]
    ) -> CommandToolBuilder:
        
        for tinput_id, history in collector.histories.items():
            if tinput_id not in valid_tinput_ids:
                if len(history.sources) == 1 and isinstance(history.sources[0], InputNodeSelector):
                    node = history.sources[0].input_node
                    if node.default is not None:
                        tinput = [x for x in tool._inputs if x.id() == tinput_id][0]  
                        tinput.default = node.default
        return tool

def get_default_tinputs(tool: CommandToolBuilder, valid_tinput_ids: set[str]) -> set[str]:
    extra_tinput_ids: set[str] = set()
    for tinput in tool._inputs:
        if tinput.id() not in valid_tinput_ids and tinput.default is not None:
            extra_tinput_ids.add(tinput.id())
    return extra_tinput_ids

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

def prune_tool(tool: CommandToolBuilder, valid_tinput_ids: set[str]) -> CommandToolBuilder:
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

def apply_pruned_tool(wf: WorkflowBuilder, pruned_tool: CommandToolBuilder) -> WorkflowBuilder:
    for step in wf.step_nodes.values():
        if isinstance(step.tool, WorkflowBuilder):
            step.tool = apply_pruned_tool(step.tool, pruned_tool)
        elif step.tool.id() == pruned_tool.id():
            step.tool = pruned_tool
        else:
            continue
    return wf