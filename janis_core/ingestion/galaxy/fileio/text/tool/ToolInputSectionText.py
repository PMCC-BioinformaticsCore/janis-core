


from dataclasses import dataclass
from typing import Tuple

from ..TextRender import TextRender
from ..tool.ToolInputText import ToolInputText
from .. import ordering

from janis_core.ingestion.galaxy.gx.command.components import InputComponent
from janis_core.ingestion.galaxy.gx.command.components import Flag
from janis_core.ingestion.galaxy.gx.command.components import Option
from janis_core.ingestion.galaxy.gx.command.components import Positional



@dataclass
class ToolInputSectionText(TextRender):
    def __init__(self, entity: list[InputComponent]):
        super().__init__()
        self.entity = entity

    @property
    def imports(self) -> list[Tuple[str, str]]:
        imports: list[Tuple[str, str]] = []
        for inp in self.entity:
            imports += ToolInputText(inp).imports
        imports = list(set(imports))
        return ordering.order_imports(imports)

    def render(self) -> str:
        out_str: str = ''

        out_str += 'inputs = ['

        out_str += '\n\t# Positionals\n'
        for positional in self._get_positionals():
            render = ToolInputText(positional)
            out_str += f'{render.render()},\n'

        out_str += '\n\t# Flags\n'
        for flag in self._get_flags():
            render = ToolInputText(flag)
            out_str += f'{render.render()},\n'

        out_str += '\n\t# Options\n'
        for option in self._get_options():
            render = ToolInputText(option)
            out_str += f'{render.render()},\n'
        
        out_str += '\n]\n'
        return out_str

    def _get_positionals(self) -> list[Positional]:
        positionals = [x for x in self.entity if isinstance(x, Positional)]
        return ordering.order_positionals(positionals)
    
    def _get_flags(self) -> list[Flag]:
        flags = [x for x in self.entity if isinstance(x, Flag)]
        return ordering.order_flags(flags)
    
    def _get_options(self) -> list[Option]:
        options = [x for x in self.entity if isinstance(x, Option)]
        return ordering.order_options(options)

    