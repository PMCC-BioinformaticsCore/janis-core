


from janis_core.ingestion.galaxy.internal_model.workflow import Workflow
from janis_core.ingestion.galaxy.internal_model.workflow.step.inputs import (
    ConnectionInputValue, 
    InputValue, 
    StaticInputValue, 
    WorkflowInputInputValue
)


def update_components_optionality(value: InputValue, janis: Workflow) -> None:
    # given a step input value which is driving a tool input (value.component is the tool input)
    # update optionality information using the source which drives this tool input
    if not value.component:
        return
    
    match value:
        case WorkflowInputInputValue():
            optionality_from_wf_input(value, janis)
        case ConnectionInputValue():
            optionality_from_connection(value, janis)
        case StaticInputValue():
            optionality_from_static_value(value)
        case _:
            raise RuntimeError

def optionality_from_wf_input(value: WorkflowInputInputValue, janis: Workflow) -> None:
    # if workflow input (src) is optional, tool input (dest) should also be optional  
    assert(value.component)
    if not value.is_runtime:
        w_inp = janis.get_input(query_uuid=value.input_uuid)
        if w_inp.optional and not value.component.optional:
            value.component.forced_optionality = True

def optionality_from_connection(value: ConnectionInputValue, janis: Workflow) -> None:
    # if prev step tool output (src) is optional, 
    # this tool input (dest) should also be optional  
    assert(value.component)
    s_out = janis.get_step_output(query_uuid=value.out_uuid) 
    if s_out.tool_output.optional and not value.component.optional:
        value.component.forced_optionality = True

def optionality_from_static_value(value: StaticInputValue) -> None:
    # if tool input is being driven by static value, 
    # just set that as the tool input default and mark as optional (now has default)
    assert(value.component)
    value.component.forced_default = value.raw_value 
    value.component.forced_optionality = True
