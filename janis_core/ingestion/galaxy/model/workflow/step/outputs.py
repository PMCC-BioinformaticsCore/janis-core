


from dataclasses import dataclass
from uuid import uuid4
from janis_core.ingestion.galaxy.gx.command.components import OutputComponent


@dataclass
class StepOutput:
    step_uuid: str
    is_wflow_out: bool
    tool_output: OutputComponent

    def __post_init__(self):
        self.uuid = str(uuid4())


