

from janis_core.ingestion.galaxy import runtime
from janis_core.ingestion.galaxy.internal_model.workflow import Workflow
from janis_core.ingestion.galaxy.internal_model.workflow import WorkflowStep
from janis_core.ingestion.galaxy.internal_model.workflow.input import WorkflowInput
from janis_core.ingestion.galaxy.internal_model.workflow.step.inputs import InputValue, WorkflowInputInputValue
from janis_core.ingestion.galaxy.datatypes.core import file_t

# from janis_core.ingestion.galaxy.gxtool import 
from janis_core.ingestion.galaxy.gxtool.model import XMLScript
from janis_core.ingestion.galaxy.gxtool.model import XMLConfigfile
from janis_core.ingestion.galaxy.gxtool.command.components import InputComponent


def handle_step_script_configfile_inputs(janis: Workflow) -> None:
    # sets tool input values as default
    # tool components which accepts scripts from the galaxy wrapper $__tool_directory
    for i_step in janis.steps:
        runtime.tool.set(from_wrapper=i_step.metadata.wrapper)
        handle_step(janis, i_step)

def handle_step(janis: Workflow, i_step: WorkflowStep) -> None:
    for script in i_step.tool.scripts:
        handle_script(janis, i_step, script)
    
    for configfile in i_step.tool.configfiles:
        handle_script(janis, i_step, configfile)

def handle_script(janis: Workflow, i_step: WorkflowStep, script: XMLScript | XMLConfigfile) -> None:
    # create worklow input & add to workflow
    winp = create_workflow_input(i_step, script)
    janis.add_input(winp)
    # create step value & add to step 
    component = get_script_component(i_step, script)
    invalue = create_workflow_value(component, winp)
    i_step.inputs.add(invalue)

def get_script_component(i_step: WorkflowStep, script: XMLScript | XMLConfigfile) -> InputComponent:
    for component in i_step.tool.inputs:
        if component.gxparam:
            if component.gxparam.name == script.varname:
                return component
    raise RuntimeError

def create_workflow_input(i_step: WorkflowStep, script: XMLScript | XMLConfigfile) -> WorkflowInput:
    return WorkflowInput(
        _name=script.varname,
        array=False,
        is_runtime=True,
        datatype=file_t,
        optional=False,
    )

def create_workflow_value(i_target: InputComponent, winp: WorkflowInput) -> InputValue:
    return WorkflowInputInputValue(
        component=i_target,
        input_uuid=winp.uuid,
        is_runtime=True
    )

# def get_wrapper_script_path(local_path: str, i_step: WorkflowStep) -> str:
#     wrapper = i_step.metadata.wrapper
#     wrapper_path = fetch_wrapper(
#         owner= wrapper.owner,
#         repo= wrapper.repo,
#         revision= wrapper.revision,
#         tool_id= wrapper.tool_id
#     )
#     wrapper_dir = wrapper_path.rsplit('/', 1)[0]
#     return f'{wrapper_dir}/{local_path}'


