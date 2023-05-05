



from typing import Any, Type
import json

from copy import deepcopy
from janis_core.ingestion.galaxy.gx.gxtool import load_xmltool
from janis_core.ingestion.galaxy import settings
from janis_core.ingestion.galaxy.gx.gxtool import XMLToolDefinition

# old
from .expand_old import expand_tool_state_old
from .flatten_old import get_flattened_tool_state_old
from .resolve_old import resolve_values_old
from .standardisation_old import standardise_tool_state_old

# from .resolve import resolve_values
# from .flatten import flatten

from .filters import (
    ReplaceNullWithVarname,
    ReplaceConnectedWithVarname,
    ReplaceRuntimeWithVarname,
    ReplaceBoolWithValue,
    IgnoreCurrentCase,
    Flatten,
    Filter
)

# MODULE ENTRY
def load_tool_state(
    step: dict[str, Any], 
    additional_filters: list[str]=[],
    ) -> dict[str, Any]:
    """ 
    loads the 'tool_state' from a galaxy step into a json dict. 
    performs optional filtering on the dict to resolve values and remove unnecessary keys.
    """
    local_filters = get_local_filters_to_apply(additional_filters)
    global_filters = get_global_filters_to_apply(additional_filters)
    xmltool = load_xmltool(settings.tool.tool_path)
    tool_state = json.loads(step['tool_state'])
    tool_state = apply_local_filters(local_filters, tool_state, xmltool)
    tool_state = apply_global_filters(global_filters, tool_state, xmltool)
    return tool_state

def apply_local_filters(filters: list[Type[Filter]], tool_state: dict[str, Any], xmltool: XMLToolDefinition) -> dict[str, Any]:
    curr_path: list[str] = []
    tool_state = do_apply_local_filters(filters, tool_state, curr_path, xmltool)
    return tool_state

def do_apply_local_filters(filters: list[Type[Filter]], tool_state: dict[str, Any], path: list[str], xmltool: XMLToolDefinition) -> dict[str, Any]:
    for fclass in filters:
        f = fclass(tool_state, path, xmltool)
        tool_state = f.apply()
    
    for key, value in tool_state.items():
        if isinstance(value, dict):
            curr_path = deepcopy(path)
            curr_path.append(key)
            tool_state[key] = do_apply_local_filters(filters, value, curr_path, xmltool)  # type: ignore
    
    return tool_state

def apply_global_filters(filters: list[Type[Filter]], tool_state: dict[str, Any], xmltool: XMLToolDefinition) -> dict[str, Any]:
    curr_path: list[str] = []
    
    for fclass in filters:
        f = fclass(tool_state, curr_path, xmltool)
        tool_state = f.apply()
    
    return tool_state

def get_local_filters_to_apply(additional_filters: list[str]) -> list[Type[Filter]]:
    all_filters = [
        ReplaceNullWithVarname,
        ReplaceConnectedWithVarname,
        ReplaceRuntimeWithVarname,
        ReplaceBoolWithValue,
        IgnoreCurrentCase,
    ]
    default_filters = [
        'ReplaceConnectedWithVarname',
        'ReplaceRuntimeWithVarname',
        'ReplaceBoolWithValue',
        'IgnoreCurrentCase',
    ]
    filter_names = set(default_filters + additional_filters)
    filters_to_apply = [x for x in all_filters if x.__name__ in filter_names]  
    return filters_to_apply  # type: ignore

def get_global_filters_to_apply(additional_filters: list[str]) -> list[Type[Filter]]:
    all_filters = [
        Flatten
    ]
    default_filters: list[str] = []
    filter_names = set(default_filters + additional_filters)
    filters_to_apply = [x for x in all_filters if x.__name__ in filter_names]  
    return filters_to_apply  # type: ignore

    


# MODULE ENTRY
def load_tool_state_dep(step: dict[str, Any]) -> dict[str, Any]:
    step['tool_state'] = expand_tool_state_old(step)  # string -> json object
    step['tool_state'] = get_flattened_tool_state_old(step)
    step['tool_state'] = resolve_values_old(step)
    step['tool_state'] = standardise_tool_state_old(step)
    return step['tool_state']

