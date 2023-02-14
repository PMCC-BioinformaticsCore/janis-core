

from typing import Any

value_translations = {
        'false': False,
        'true': True,
    }


def standardise_tool_state(gxstep: dict[str, Any]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for key, val in gxstep['tool_state'].items():
        out[key] = standardise_tool_state_value(val)
    return out

def standardise_tool_state_value(value: Any) -> Any:
    value = handle_array(value)
    value = handle_translations(value)
    return value

def handle_array(value: Any) -> Any:
    if isinstance(value, list):
        return None
        # if len(value) == 0:
        #     value = None
        # else:
        #     value = value[0]
    return value

def handle_translations(value: Any) -> Any:
    if value in value_translations:
        return value_translations[value]
    return value
