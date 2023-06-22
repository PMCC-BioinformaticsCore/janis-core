


from typing import Any


from galaxy.tool_util.verify.asserts.text import (
    assert_has_line,
    assert_has_line_matching,
    assert_has_text,
    assert_has_text_matching,
    assert_not_has_text
)
from galaxy.tool_util.verify.asserts.tabular import (
    assert_has_n_columns
)


"""
# Implemented
'lines_diff': lines_diff_strategy,
'has_size': has_size_strategy,
'has_text': has_text_strategy,
'has_n_lines': not_implemented_strategy,
'has_n_columns': not_implemented_strategy,

# not Implemented
'ftype': ftype_strategy,
'not_has_text': not_implemented_strategy,
'has_text_matching': not_implemented_strategy,
'has_line': not_implemented_strategy,
'has_line_matching': not_implemented_strategy,
"""


def has_text(output_value: Any, kwargs: dict[str, Any]) -> bool:
    try:
        assert_has_text(output_value, kwargs['expected_value'], n=None)
    except AssertionError:
        return False
    return True

def has_n_columns(output_value: Any, kwargs: dict[str, Any]) -> bool:
    try:
        assert_has_n_columns(output_value, kwargs['expected_value'])
    except AssertionError:
        return False
    return True


def has_text_matching(output_value: Any, kwargs: dict[str, Any]) -> bool:
    try:
        assert_has_text_matching(output_value, kwargs['expected_value'])
    except AssertionError:
        return False
    return True


