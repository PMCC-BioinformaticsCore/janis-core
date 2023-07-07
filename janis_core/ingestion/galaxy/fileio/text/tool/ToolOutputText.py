

from dataclasses import dataclass
from typing import Optional, Tuple

from janis_core.ingestion.galaxy.gxtool.command.components import OutputComponent
from janis_core.ingestion.galaxy.gxtool.command.components import InputOutput
from janis_core.ingestion.galaxy.gxtool.command.components import RedirectOutput
from janis_core.ingestion.galaxy.gxtool.command.components import WildcardOutput

from ..TextRender import TextRender
from .. import formatting
from .. import ordering

from janis_core.ingestion.galaxy import tags  # remove if possible


def format_selector_str(output: OutputComponent) -> Optional[str]:
    match output: 
        case RedirectOutput():
            return None
        case InputOutput():
            input_comp_uuid = output.input_component.uuid
            input_comp_tag = tags.get(input_comp_uuid)
            return f'InputSelector("{input_comp_tag}")'
        case WildcardOutput():
            pattern = 'unknown'
            if hasattr(output.gxparam, 'from_work_dir') and output.gxparam.from_work_dir:
                pattern = output.gxparam.from_work_dir
            elif output.gxparam.discover_pattern:
                pattern = output.gxparam.discover_pattern
            return f'WildcardSelector(r"{pattern}")'
        case _:
            pass


@dataclass
class ToolOutputText(TextRender):
    def __init__(self, entity: OutputComponent):
        super().__init__()
        self.entity = entity

    @property
    def imports(self) -> list[Tuple[str, str]]:
        jtype = self.entity.datatype
        imports: list[Tuple[str, str]] = []
        imports.append((jtype.import_path, jtype.classname))

        # Stdout class import
        if isinstance(self.entity, RedirectOutput):
            imports.append(('janis_core', 'Stdout'))

        # Selector class import
        selector_str = format_selector_str(self.entity)
        if selector_str:
            selector = selector_str.split('(', 1)[0]
            imports.append(('janis_core', selector))
        
        # Array class import
        if self.entity.array:
            imports.append(('janis_core', 'Array'))

        imports = list(set(imports))
        return ordering.order_imports(imports)

    def render(self) -> str:
        e = self.entity
        selector_str = format_selector_str(e)
        doc = formatting.format_docstring(e)
        datatype_str = formatting.format_typestr(e)
        if isinstance(e, RedirectOutput):
            datatype_str = f'Stdout'
        out_str: str = ''
        out_str += '\tToolOutput(\n'
        out_str += f"\t\t'{tags.get(e.uuid)}',\n"
        out_str += f"\t\t{datatype_str},\n"
        if selector_str:
            out_str += f"\t\tselector={selector_str},\n" 
        out_str += f'\t\tdoc="{doc}",\n' if doc else ''
        out_str += '\t)'
        return out_str

    