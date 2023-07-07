


from typing import Any


def sanitise_default_value(value: Any) -> Any:
    value = _sanitise_env_var(value)
    return value

def _sanitise_env_var(value: Any) -> Any:
    if isinstance(value, str) and value.startswith('$'):
        value = None
    return value