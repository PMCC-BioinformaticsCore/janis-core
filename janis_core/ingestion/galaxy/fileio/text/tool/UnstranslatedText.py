

from __future__ import annotations
from typing import TYPE_CHECKING, Tuple

if TYPE_CHECKING:
    from janis_core.ingestion.galaxy.model.workflow import WorkflowStep

from ..TextRender import TextRender


class UntranslatedText(TextRender):
    def __init__(self, entity: WorkflowStep):
        super().__init__()
        self.entity = entity

    @property
    def imports(self) -> list[Tuple[str, str]]:
        return []

    def render(self) -> str:
        step = self.entity
        out_str: str = ''
        if step.preprocessing:
            out_str += '\n# PRE-PROCESSING ---------------------\n'
            out_str += f'{step.preprocessing}\n'
        if step.postprocessing:
            out_str += '\n# POST-PROCESSING ---------------------\n'
            out_str += f'{step.postprocessing}\n'
        return out_str