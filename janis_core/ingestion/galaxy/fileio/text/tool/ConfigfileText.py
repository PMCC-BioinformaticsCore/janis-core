

from typing import Tuple
from janis_core.ingestion.galaxy.gxtool.model import XMLConfigfile

from ..TextRender import TextRender


class ConfigfileText(TextRender):
    def __init__(self, entity: XMLConfigfile):
        super().__init__()
        self.entity = entity

    @property
    def imports(self) -> list[Tuple[str, str]]:
        return []

    def render(self) -> str:
        return self.entity.contents
    
