

from typing import Any, Optional, Protocol, Tuple

from janis_core.ingestion.galaxy.gx.command.components import CommandComponent
from janis_core.ingestion.galaxy.gx.command.components import InputComponent
from janis_core.ingestion.galaxy.gx.command.components import Option

from janis_core.ingestion.galaxy.datatypes import JanisDatatype



def format_docstring(component: CommandComponent) -> Optional[str]:
    raw_doc = component.docstring
    if raw_doc:
        return raw_doc.replace('"', "'")
    return None

def get_wrapped_default(component: InputComponent) -> Optional[str]:
    default = component.default_value
    if isinstance(component, Option) and default is not None:
        if isinstance(default, str) and '$' in default:
            default = None
    if should_quote(default, component):
        return f'"{default}"'
    return default

def format_imports(imports: list[Tuple[str, str]]) -> str:
    out_str: str = ''
    for x in imports:
        out_str += f'from {x[0]} import {x[1]}\n'
    return out_str


class TypeFormattable(Protocol):
    @property
    def optional(self) -> bool:
        ...

    @property
    def array(self) -> bool:
        ...
    
    @property
    def datatype(self) -> JanisDatatype:
        ...
    

def should_quote(default: Any, component: TypeFormattable) -> bool:
    jtype = component.datatype
    if jtype.classname == 'Int' or jtype.classname == 'Float':
        return False
    elif default is None:
        return False
    return True

def datatype_permits_default(jtype: JanisDatatype) -> bool:
    # check this component's datatypes permit a default value
    types_with_allowed_default = ['String', 'Float', 'Int', 'Boolean', 'Double']
    if jtype.classname in types_with_allowed_default:
        return True
    return False

def format_typestr(entity: TypeFormattable, fmt: str='definition') -> str:
    if fmt == 'definition':
        return fmt_typestr_definition(entity)
    elif fmt == 'value':
        return fmt_typestr_value(entity)
    else:
        raise RuntimeError()

def fmt_typestr_definition(entity: TypeFormattable) -> str:
    jtype = entity.datatype
    # not array not optional
    if not entity.optional and not entity.array: 
        out_str = f'{jtype.classname}'
    
    # array and not optional
    elif not entity.optional and entity.array: 
        out_str = f'Array({jtype.classname})'
    
    # not array and optional
    elif entity.optional and not entity.array: 
        out_str = f'{jtype.classname}(optional=True)'
    
    # array and optional
    elif entity.optional and entity.array: 
        out_str = f'Array({jtype.classname}(), optional=True)'
    return out_str

    # if len(types) > 1:
    #     dtype = ', '.join([x.classname for x in types])
    #     dtype = "UnionType(" + dtype + ")"
    # else:
    #     dtype = types[0].classname

def fmt_typestr_value(entity: TypeFormattable) -> str:
    jtype = entity.datatype
    # not array not optional
    if not entity.optional and not entity.array: 
        out_str = f'{jtype.classname}'
    
    # array and not optional
    elif not entity.optional and entity.array: 
        out_str = f'{jtype.classname} ARRAY'
    
    # not array and optional
    elif entity.optional and not entity.array: 
        out_str = f'{jtype.classname} OPTIONAL'
    
    # array and optional
    elif entity.optional and entity.array: 
        out_str = f'{jtype.classname} ARRAY OPTIONAL'
    return out_str

    # if len(types) > 1:
    #     dtype = ', '.join([x.classname for x in types])
    #     dtype = "UnionType(" + dtype + ")"
    # else:
    #     dtype = types[0].classname

    

