

from janis_core.ingestion.galaxy import settings
from janis_core.ingestion.galaxy import paths
import shutil
from janis_core.ingestion.galaxy import expressions

from janis_core.ingestion.galaxy.model.workflow import Workflow
from janis_core.ingestion.galaxy.model.workflow import WorkflowStep
from janis_core.ingestion.galaxy.model.workflow.input import WorkflowInput
from janis_core.ingestion.galaxy.model.workflow.step.inputs import InputValue, WorkflowInputInputValue
from janis_core.ingestion.galaxy.datatypes.core import file_t

from ...command.components import InputComponent
from ...command.components import Positional
from ...command.components import Option
from ...wrappers import fetch_wrapper


def handle_step_script_inputs(janis: Workflow) -> None:
    # sets tool input values as default
    # tool components which accepts scripts from the galaxy wrapper $__tool_directory
    for j_step in janis.steps:
        settings.tool.set(from_wrapper=j_step.metadata.wrapper)
        handle_step(janis, j_step)

def handle_step(janis: Workflow, j_step: WorkflowStep) -> None:
    for component in get_linkable_components(j_step):
        if component_is_script(component):
            # get filepath for script
            galaxy_path: str = component.values.unique[0] # type: ignore
            local_path = galaxy_path.replace('$__tool_directory__/', '')
            dest = paths.script(local_path)
            src = get_wrapper_script_path(local_path, j_step)

            # copy the script to the destination folder
            shutil.copy2(src, dest)
            
            # create worklow input & add to workflow
            winp = create_workflow_input(j_step, component, dest)
            janis.add_input(winp)
            
            # create step value & add to step 
            invalue = create_workflow_value(component, winp)
            j_step.inputs.add(invalue)

def component_is_script(component: InputComponent) -> bool:
    if isinstance(component, Positional) or isinstance(component, Option):
        if len(component.values.unique) == 1:
            if expressions.is_script(component.values.unique[0]):
                return True
    return False

def get_wrapper_script_path(local_path: str, step: WorkflowStep) -> str:
    wrapper = step.metadata.wrapper
    wrapper_path = fetch_wrapper(
        owner= wrapper.owner,
        repo= wrapper.repo,
        revision= wrapper.revision,
        tool_id= wrapper.tool_id
    )
    wrapper_dir = wrapper_path.rsplit('/', 1)[0]
    return f'{wrapper_dir}/{local_path}'

def get_linkable_components(j_step: WorkflowStep) -> list[InputComponent]:
    out: list[InputComponent] = []
    for component in j_step.tool.inputs:
        if not j_step.inputs.get(component.uuid):
            out.append(component)
    return out

def create_workflow_input(j_step: WorkflowStep, j_target: InputComponent, filepath: str) -> WorkflowInput:
    step_tag = j_step.tag
    input_tag = j_target.tag
    return WorkflowInput(
        _name=f'{step_tag}_{input_tag}',
        array=j_target.array,
        is_runtime=True,
        datatype=file_t,
        optional=j_target.optional,
        value=filepath
    )

def create_workflow_value(j_target: InputComponent, winp: WorkflowInput) -> InputValue:
    return WorkflowInputInputValue(
        component=j_target,
        input_uuid=winp.uuid,
        is_runtime=True
    )
