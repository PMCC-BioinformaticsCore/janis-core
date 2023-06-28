


from dataclasses import dataclass
from typing import Tuple
from janis_core.ingestion.galaxy.gx.command.components import OutputComponent

from ..TextRender import TextRender
from ..tool.ToolOutputText import ToolOutputText
from .. import ordering


@dataclass
class ToolOutputSectionText(TextRender):
    def __init__(self, entity: list[OutputComponent]):
        super().__init__()
        self.entity = entity
    
    @property
    def imports(self) -> list[Tuple[str, str]]:
        imports: list[Tuple[str, str]] = []
        for out in self.entity:
            imports += ToolOutputText(out).imports
        imports = list(set(imports))
        return ordering.order_imports(imports)

    def render(self) -> str:
        out_str: str = ''
        
        out_str += 'outputs = [\n'
        for out in self.entity:
            render = ToolOutputText(out)
            out_str += f'\t{render.render()},\n'
        out_str += '\n]\n'

        return out_str