

from janis_core import WorkflowBuilder, CommandToolBuilder
from janis_core.workflow.workflow import StepNode
from janis_core.operators.selectors import InputSelector
from janis_core.translations.common import trace
from .history import TaskInputCollector
from janis_core import translation_utils as utils
from janis_core.translation_utils import DTypeType

def prune_tools_and_sources(main_wf: WorkflowBuilder, tools: dict[str, CommandToolBuilder]) -> None:
    for tool in tools.values():
        prune(main_wf, tool)

def prune(main_wf: WorkflowBuilder, tool: CommandToolBuilder) -> None:
    # for each tool input, get all sources from each wf step which feed this input 
    collector = TaskInputCollector(tool)
    collector.collect(main_wf)

    ### STEP INPUTS ###
    # identifying which tinputs need to be fed via step.sources
    valid_tinput_ids: set[str] = set() 
    # identify tinputs which need to be kept based on step.sources
    valid_tinput_ids = valid_tinput_ids | get_step_referenced_tinputs(collector)
    # identify tinputs which optional files with default referencing another input
    valid_tinput_ids = valid_tinput_ids | get_optional_file_inputs_w_default(tool)
    # identify tinputs which are referenced in tool outputs
    valid_tinput_ids = valid_tinput_ids | get_output_referenced_tinputs(main_wf, tool)

    ### STATIC VALUES & SOURCES ###
    # identify tinputs which have a single static value as source & migrate the source to tool default
    migrate_single_statics_to_defaults(valid_tinput_ids, collector, tool)
    # remove step.sources which are not needed, migrate inputnode srcs to raw values
    prune_sources(main_wf, tool.id(), valid_tinput_ids)

    # add tinputs which reference a previously validated tinput
    valid_tinput_ids = valid_tinput_ids | get_tinput_reference_tinputs(tool, valid_tinput_ids)
    # add tinputs which have a non-null default
    valid_tinput_ids = valid_tinput_ids | get_default_tinputs(tool, valid_tinput_ids)
    
    # remove the ToolInputs / ToolArguments which are not needed
    prune_tool(main_wf, tool.id(), valid_tinput_ids)
    # apply the pruned tool to the workflow
    # apply_pruned_tool(main_wf, tool)

def get_step_referenced_tinputs(collector: TaskInputCollector) -> set[str]:
    # get the tinputs which are needed based on step inputs
    filtered_tinputs: set[str] = set()
    
    for tinput_id, history in collector.histories.items():
        # RULE 1: mandatory tool inputs are always kept
        if not history.is_optional:  
            filtered_tinputs.add(tinput_id)
        # RULE 2: tool inputs which are fed by mandatory types are kept
        elif len(history.mandatory_input_sources) >= 1:
            filtered_tinputs.add(tinput_id)
        # RULE 3: if has multiple sources, keep
        elif len(history.sources) >= 2:
            filtered_tinputs.add(tinput_id)
        # RULE 4: if has step connections, keep
        elif len(history.connection_sources) >= 1:
            filtered_tinputs.add(tinput_id)
        # RULE 5: if weird sources, keep
        elif len(history.other_sources) >= 1:
            filtered_tinputs.add(tinput_id)
        # RULE 6: (edge case) if tool used in 2+ steps, but tinput has only 1 source, keep.
        #         this is needed because if we are driving the tinput's value from a source in 
        #         one step, then it's value is driven by its default in another. they could be different. 
        elif len(history.sources) == 1 and collector.step_count >= 2:
            filtered_tinputs.add(tinput_id)
        else:
            continue
        # # RULE 5: if only has single placeholder source, ignore
        # elif len(history.sources) == 1 and len(history.placeholder_sources) == 1:
        #     continue
        # # RULE 6: if has 1+ sources which are workflow inputs, keep
        # elif len(history.input_sources) >= 1:
        #     tinputs_to_keep.add(tinput_id)
        # else:
        #     continue
    
    return filtered_tinputs

def get_optional_file_inputs_w_default(tool: CommandToolBuilder) -> set[str]:
    filtered_tinputs: set[str] = set()
    for tinput in tool._inputs:
        dtt = utils.get_dtt(tinput.input_type) # type: ignore
        if dtt == DTypeType.FILE:
            if tinput.input_type.optional and tinput.default is not None:
                filtered_tinputs.add(tinput.id())
    return filtered_tinputs

def get_output_referenced_tinputs(main_wf: WorkflowBuilder, tool: CommandToolBuilder) -> set[str]:
    # TODO check if the output is consumed as a source somewhere else in the workflow. 
    # if not we can remove it. 
    filtered_tinputs: set[str] = set()
    for tout in tool._outputs:
        ref_vars = trace.trace_referenced_variables(tout, tool)
        filtered_tinputs = filtered_tinputs | ref_vars
    return filtered_tinputs

def migrate_single_statics_to_defaults(
    valid_tinput_ids: set[str], 
    collector: TaskInputCollector, 
    tool: CommandToolBuilder
    ) -> None:
    for tinput_id, history in collector.histories.items():
        if tinput_id in valid_tinput_ids:
            continue
        if len(history.sources) == 1 and len(history.placeholder_sources) == 1:
            node = history.input_sources[0].value.input_node
            tinput = [x for x in tool._inputs if x.id() == tinput_id][0] 
            tinput.default = node.default

def prune_sources(local_wf: WorkflowBuilder, tool_id: str, valid_tinput_ids: set[str]) -> None:
    for step in local_wf.step_nodes.values():
        if isinstance(step.tool, WorkflowBuilder):
            prune_sources(step.tool, tool_id, valid_tinput_ids)
        elif isinstance(step.tool, CommandToolBuilder) and step.tool.id() == tool_id:
            do_prune_sources(step, valid_tinput_ids)
        else:
            continue

def do_prune_sources(step: StepNode, valid_tinput_ids: set[str]) -> None:
    # remove sources which are not needed
    invalid_sources = set(step.sources.keys()) - valid_tinput_ids
    for tinput_id in invalid_sources:
        del step.sources[tinput_id]

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

def prune_tool(local_wf: WorkflowBuilder, tool_id: str, valid_tinput_ids: set[str]) -> None:
    for step in local_wf.step_nodes.values():
        if isinstance(step.tool, WorkflowBuilder):
            prune_tool(step.tool, tool_id, valid_tinput_ids)
        elif isinstance(step.tool, CommandToolBuilder) and step.tool.id() == tool_id:
            do_prune_tool_inputs(step.tool, valid_tinput_ids)
            do_prune_tool_arguments(step.tool, valid_tinput_ids)
        else:
            continue
    
def do_prune_tool_inputs(tool: CommandToolBuilder, valid_tinput_ids: set[str]) -> None:
    # early exit
    if not tool._inputs:
        return 
    
    items_to_delete: list[int] = []
    for i, tinput in enumerate(tool._inputs):     # type: ignore
        if tinput.id() not in valid_tinput_ids:
            items_to_delete.append(i)
    
    for i in sorted(items_to_delete, reverse=True):
        del tool._inputs[i]     # type: ignore

def do_prune_tool_arguments(tool: CommandToolBuilder, valid_tinput_ids: set[str]) -> None:
    # early exit
    if not tool._arguments:
        return 
    
    items_to_delete: list[int] = []
    for i, targ in enumerate(tool._arguments):     # type: ignore
        refs = trace.trace_entities(targ, tool)
        input_refs = [ref for ref in refs if isinstance(ref, InputSelector)]
        if not input_refs:
            continue

        if not all([ref.input_to_select in valid_tinput_ids for ref in input_refs]):
            items_to_delete.append(i)
    
    for i in sorted(items_to_delete, reverse=True):
        del tool._arguments[i]     # type: ignore
    

