
from typing import Tuple

from janis_core.ingestion.galaxy.internal_model.workflow import Workflow
from ..TextRender import TextRender


class ConfigText(TextRender):
    def __init__(self, entity: Workflow):
        super().__init__()
        self.entity = entity

    @property
    def imports(self) -> list[Tuple[str, str]]:
        raise NotImplementedError()

    def render(self) -> str:
        raise NotImplementedError()
    
