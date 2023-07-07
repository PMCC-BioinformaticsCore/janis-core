




from janis_core.ingestion.galaxy.internal_model.workflow import Workflow
from janis_core.ingestion.galaxy.internal_model.workflow.step.inputs import InputValue


def update_components_array(value: InputValue, janis: Workflow) -> None:
    """
    correct information about tool components using InputValue.
    ie if workflow input is array, but tool component is not, 
    update tool component to be an array. 
    """
    pass
    # if value.component:  # if linked to tool component
    #     if isinstance(value, WorkflowInputInputValue):
    #         if not value.is_runtime:
    #             w_inp = janis.get_input(query_uuid=value.input_uuid)
    #             value.component.forced_array = w_inp.array
    #     elif isinstance(value, ConnectionInputValue):
    #         s_out = janis.get_step_output(query_uuid=value.out_uuid) # type: ignore
    #         value.component.forced_array = s_out.tool_output.array

