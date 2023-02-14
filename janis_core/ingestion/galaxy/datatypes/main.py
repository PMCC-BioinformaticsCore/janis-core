

from typing import Optional, Any

from .JanisDatatype import JanisDatatype
from .strategies import strategy_map
from .conversion import janis_to_core
from .conversion import select_primary_core_type
from .register import register


def get(entity: Any, entity_type: Optional[str]=None, source: str='galaxy') -> JanisDatatype:
    etype = entity_type if entity_type else entity.__class__.__name__
    strategy = strategy_map[etype] 
    jtypes = strategy().get(entity)
    core_types = janis_to_core(jtypes)
    return select_primary_core_type(core_types)

def populate() -> None:
    register.populate()

