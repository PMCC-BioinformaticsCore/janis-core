

from typing import Callable
from .filters import (
    flatten_nesting,
    flatten_multiline_strings,
    translate_variable_markers,
    standardise_variable_format,
    simplify_sh_constructs,
    simplify_galaxy_dynamic_vars,
    remove_cheetah_comments,
    replace_function_calls,
    replace_backticks,
    remove_empty_quotes,
    interpret_raw
)



def simplify_cmd(text: str, purpose: str) -> str:
    simplifier_map = {
        'templating': CheetahSimplifier(),
        'parsing': ParsingSimplifier(),
    }

    simplifier = simplifier_map[purpose]
    return simplifier.simplify(text)


# classes 

class CommandSimplifier:
    filters: list[Callable[[str], str]] = []

    def simplify(self, cmdstr: str) -> str:
        for filter_func in self.filters:
            cmdstr = filter_func(cmdstr)
        return cmdstr

class CheetahSimplifier(CommandSimplifier):
    filters: list[Callable[[str], str]] = [
        simplify_galaxy_dynamic_vars,
    ]

# class TestSimplifier(CommandSimplifier):
#     filters: list[Callable[[str], str]] = [
#         translate_variable_markers,
#         standardise_variable_format,
#         simplify_galaxy_dynamic_vars,
#         simplify_sh_constructs,
#         replace_backticks,
#         remove_empty_quotes,
#     ]

class ParsingSimplifier(CommandSimplifier):
    filters: list[Callable[[str], str]] = [
        flatten_nesting,
        remove_cheetah_comments,
        flatten_multiline_strings,
        replace_function_calls,
        replace_backticks,
        standardise_variable_format,  # ?
        simplify_sh_constructs,
        simplify_galaxy_dynamic_vars,
        # remove_empty_quotes,
        # interpret_raw # ?
    ]

