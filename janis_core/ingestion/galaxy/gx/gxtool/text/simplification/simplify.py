

from typing import Callable
from .filters import (
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



def simplify_cmd(text: str, output: str) -> str:
    simplifier_map = {
        'test': TestSimplifier(),
        'xml': XMLSimplifier(),
        'cheetah': PartialCheetahEvalSimplifier(),
    }
    simplifier = simplifier_map[output]
    return simplifier.simplify(text)



# classes 

class CommandSimplifier:
    filters: list[Callable[[str], str]] = []

    def simplify(self, cmdstr: str) -> str:
        return self.map_filters(cmdstr)

    def map_filters(self, cmdstr: str) -> str:
        for filter_func in self.filters:
            cmdstr = filter_func(cmdstr)
        return cmdstr


class PartialCheetahEvalSimplifier(CommandSimplifier):
    filters: list[Callable[[str], str]] = [
        remove_cheetah_comments,
        simplify_galaxy_dynamic_vars,
    ]


class TestSimplifier(CommandSimplifier):
    filters: list[Callable[[str], str]] = [
        translate_variable_markers,
        standardise_variable_format,
        simplify_galaxy_dynamic_vars,
        simplify_sh_constructs,
        replace_backticks,
        remove_empty_quotes,
    ]


class XMLSimplifier(CommandSimplifier):
    filters: list[Callable[[str], str]] = [
        remove_cheetah_comments,
        flatten_multiline_strings,
        replace_function_calls,
        replace_backticks,
        standardise_variable_format,  # ?
        simplify_sh_constructs,
        simplify_galaxy_dynamic_vars,
        remove_empty_quotes,
        interpret_raw
    ]

