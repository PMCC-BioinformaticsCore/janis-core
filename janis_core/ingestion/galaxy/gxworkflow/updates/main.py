

from janis_core.ingestion.galaxy.internal_model.workflow.workflow import Workflow

from .optionality import update_components_optionality
from .arrays import update_components_array
from .datatypes import update_components_datatype


def update_component_knowledge(janis: Workflow) -> None:
    for step in janis.steps:
        for value in step.inputs.linked:
            update_components_optionality(value, janis)
            update_components_array(value, janis)
            update_components_datatype(value, janis)
        