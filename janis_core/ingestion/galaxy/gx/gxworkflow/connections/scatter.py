



from janis_core.ingestion.galaxy.model.workflow.input import WorkflowInput
from janis_core.ingestion.galaxy.model.workflow.step.step import WorkflowStep
from janis_core.ingestion.galaxy.model.workflow.workflow import Workflow


def handle_scattering(janis: Workflow) -> None:
    """
    Identify the inputs in each step which need to be scattered.
    Workflow inputs which are arrays trigger any downstream steps
    to be scattered. 
    Still unsure if array outputs should also trigger scattering.
    """
    # get workflow inputs which are arrays
    for winp in janis.inputs:
        if winp.array: 
            # scatter direct child steps, then recursively scatter
            # child steps which connect to affected steps
            children = janis.get_input_children(winp.uuid)
            for child in children:
                update_scatter_inputs(child, parent=winp)
                recursive_scatter(child, janis)

def recursive_scatter(step: WorkflowStep, janis: Workflow) -> None:
    """
    Identify inputs which will be scattered on.
    Then do for child steps.
    """
    children = janis.get_step_children(step.uuid)
    for child in children:
        update_scatter_inputs(child, parent=step)
        recursive_scatter(child, janis)

def update_scatter_inputs(step: WorkflowStep, parent: WorkflowInput | WorkflowStep):
    match parent:
        case WorkflowInput():
            invalues = step.inputs.workflow_inputs
            for val in invalues:
                if val.input_uuid == parent.uuid:
                    val.scatter = True
        case WorkflowStep():
            invalues = step.inputs.connections
            for val in invalues:
                if val.step_uuid == parent.uuid:
                    val.scatter = True
        case _:
            raise RuntimeError()

