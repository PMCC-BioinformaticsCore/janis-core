
import regex as re

from .patterns import (
    INTEGER,
    FLOAT,
    SCRIPT,
    VARIABLES_FMT1,
    VARIABLES_FMT2,
)

from .matches import get_matches

def is_int(the_string: str) -> bool:
    matches = get_matches(the_string, INTEGER)
    if matches:
        return True
    return False

def is_float(the_string: str) -> bool:
    matches = get_matches(the_string, FLOAT)
    if matches:
        return True
    return False

def is_var(the_string: str) -> bool:
    matches1 = get_matches(the_string, VARIABLES_FMT1)
    matches2 = get_matches(the_string, VARIABLES_FMT2)
    if matches1 or matches2:
        return True
    return False

def is_script(the_string: str) -> bool:
    matches = get_matches(the_string, SCRIPT)
    if matches:
        return True
    return False

def has_var(the_string: str) -> bool:
    if '$' in the_string:
        return True
    return False

def is_present(word: str, text: str) -> bool:
    pattern = rf'\s{word}\s'
    if re.findall(pattern, text):
        return True
    return False
