


from abc import ABC, abstractmethod
from typing import Any
from janis_core.ingestion.galaxy import expressions
import json

# from janis_core.ingestion.galaxy.gx.command.components.inputs.Option import Option

from .JanisDatatype import JanisDatatype
from .conversion import galaxy_to_janis
from . import core



### SHARED FUNCTIONS ###

def types_from_param(entity: Any) -> list[str]:
    types = []
    if entity.gxparam and entity.gxparam.formats:
        types = entity.gxparam.formats
    
    # edge case: Option with boolean gxparam value
    # eg ASSUME_SORTED=${assume_sorted}, where ${assume_sorted} is bool param
    if entity.__class__.__name__ == 'Option' and types == ['boolean']:
        return []
    else:
        return types

def types_from_values(values: list[Any]) -> list[str]:
    if all([expressions.is_script(v) for v in values]):
        raise RuntimeError
    elif all([expressions.is_int(v) for v in values]):
        return ['integer']
    elif all([expressions.is_float(v) for v in values]):
        return ['float']
    elif all([expressions.is_var(v) for v in values]):
        return ['string']
    return ['string']

def types_from_default(entity: Any) -> list[str]:
    default = str(entity.default_value)
    if default:
        if expressions.is_int(default):
            return ['integer']
        if expressions.is_float(default):
            return ['float']
        if expressions.is_var(default) or expressions.has_var(default):
            return ['string']
    return []

def types_from_extension(output: Any) -> list[str]:
    if output.gxparam and hasattr(output.gxparam, 'discover_pattern'):
        pattern: str = output.gxparam.discover_pattern
        if pattern and '.' in pattern:
            compressed_types = ['gz', 'bz2']
            p_split = pattern.split('.')
            while p_split[-1] in compressed_types:
                p_split = p_split[:-1]
            ext = p_split[-1]
            return [f'.{ext}']
    return []


### STRATEGIES ###

class DatatypeGetStrategy(ABC):
    @abstractmethod
    def get(self, entity: Any) -> list[JanisDatatype]:
        """parses an entity to return the janis datatype(s)"""
        ...

class PositionalStrategy(DatatypeGetStrategy):
    def get(self, entity: Any) -> list[JanisDatatype]:
        gxtypes = types_from_param(entity)
        if not gxtypes:
            gxtypes = types_from_default(entity)
        if not gxtypes:
            gxtypes = types_from_values(entity.values.unique)
        return galaxy_to_janis(gxtypes)

class FlagStrategy(DatatypeGetStrategy):
    def get(self, entity: Any) -> list[JanisDatatype]:
        return [core.bool_t]

class OptionStrategy(DatatypeGetStrategy):
    def get(self, entity: Any) -> list[JanisDatatype]:
        gxtypes = types_from_param(entity)
        if not gxtypes:
            gxtypes = types_from_default(entity)
        if not gxtypes:
            gxtypes = types_from_values(entity.values.unique)
        return galaxy_to_janis(gxtypes)

class OutputStrategy(DatatypeGetStrategy):
    def get(self, entity: Any) -> list[JanisDatatype]:
        gxtypes = types_from_param(entity)
        if not gxtypes:
            gxtypes = types_from_extension(entity)
        if not gxtypes:
            gxtypes = ['file']
        return galaxy_to_janis(gxtypes)

class WorkflowInputStrategy(DatatypeGetStrategy):
    def get(self, entity: Any) -> list[JanisDatatype]:
        return [entity.datatype]

class GalaxyInputStepStrategy(DatatypeGetStrategy):
    def get(self, entity: dict[str, Any]) -> list[JanisDatatype]:
        tool_state = json.loads(entity['tool_state'])
        if 'format' in tool_state:
            return galaxy_to_janis(tool_state['format'])
        else:
            return [core.file_t]



# these are the entities which we can get a datatype for

strategy_map = {
    'Positional': PositionalStrategy,
    'Flag': FlagStrategy,
    'Option': OptionStrategy,
    'RedirectOutput': OutputStrategy,
    'InputOutput': OutputStrategy,
    'WildcardOutput': OutputStrategy,
    'GalaxyInputStep': GalaxyInputStepStrategy,
    'WorkflowInput': WorkflowInputStrategy,
}





        
