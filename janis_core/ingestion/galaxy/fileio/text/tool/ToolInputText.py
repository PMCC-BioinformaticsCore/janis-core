


from typing import Optional, Tuple
from dataclasses import dataclass

from janis_core.ingestion.galaxy.gxtool.command.components import Flag
from janis_core.ingestion.galaxy.gxtool.command.components import InputComponent
from janis_core.ingestion.galaxy.gxtool.command.components import Positional
from janis_core.ingestion.galaxy.gxtool.command.components import Option

from ..TextRender import TextRender
from .. import formatting
from .. import ordering


@dataclass
class ToolInputText(TextRender):
    def __init__(self, entity: InputComponent):
        super().__init__()
        self.entity = entity

    @property
    def imports(self) -> list[Tuple[str, str]]:
        jtype = self.entity.datatype
        imports: list[Tuple[str, str]] = []
        imports.append((jtype.import_path, jtype.classname))

        if self.entity.array:
            imports.append(('janis_core', 'Array'))
        
        imports = list(set(imports))
        return ordering.order_imports(imports)

    def render(self) -> str:
        e = self.entity
        prefix = self.format_component_prefix()
        kv_space = self.format_separation()
        doc = formatting.format_docstring(e)
        datatype_str = formatting.format_typestr(e)
        out_str: str = ''
        out_str += '\tToolInput(\n'
        out_str += f"\t\t'{e.tag}',\n"
        out_str += f"\t\t{datatype_str},\n" 
        out_str += f"\t\tprefix='{prefix}',\n" if prefix else ''
        out_str += f"\t\tseparate_value_from_prefix={kv_space},\n" if kv_space == False else ''
        out_str += f"\t\tposition={e.cmd_pos},\n"
        out_str += f"\t\tdefault={formatting.get_wrapped_default(e)},\n"
        out_str += f'\t\tdoc="{doc}",\n' if doc else ''
        out_str += '\t)'
        return out_str

    def format_component_prefix(self) -> Optional[str]:
        e = self.entity
        match e:
            case Positional():
                return None
            case Flag():
                return e.prefix
            case Option():
                return e.prefix if e.delim == ' ' else e.prefix + e.delim
            case _:
                raise RuntimeError
    
    def format_separation(self) -> bool:
        if isinstance(self.entity, Option):
            if self.entity.delim == ' ':
                return True
            else:
                return False
        else:
            return True


