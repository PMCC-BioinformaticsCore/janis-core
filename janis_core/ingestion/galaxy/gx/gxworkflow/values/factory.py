
from typing import Any, Optional

from janis_core.ingestion.galaxy.gx.command.components import InputComponent

from janis_core.ingestion.galaxy.model.workflow import (
    ConnectionInputValue, 
    WorkflowInputInputValue,
    StaticInputValue
)


def static(component: Optional[InputComponent], value: Any, default: bool=False) -> StaticInputValue:
    value = str(value)      # yeaaaaaa um yea. dw about this.
    return StaticInputValue(
        component=component,
        str_value=value,
        is_default=default
    )

def connection(component: Optional[InputComponent], step_uuid: str, out_uuid: str) -> ConnectionInputValue:
    return ConnectionInputValue(
        component=component,
        step_uuid=step_uuid,
        out_uuid=out_uuid
    )

def workflow_input(component: Optional[InputComponent], input_uuid: str, is_runtime: bool=False) -> WorkflowInputInputValue:
    return WorkflowInputInputValue(
        component=component,
        input_uuid=input_uuid,
        is_runtime=is_runtime
    )


