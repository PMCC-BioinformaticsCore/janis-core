

from typing import Any

from ..parsing.main import load_xmltool
from ..text.simplification.aliases import resolve_aliases
from ..text.simplification.main_statement import mark_main_statement
from ..text.cheetah.evaluation import sectional_evaluate
from ..text.simplification.simplify import simplify_cmd

from janis_core.ingestion.galaxy import runtime


def load_vanilla_command_str() -> str:
    """
    loads <command> section of the active tool's XML for analysis.
    simplifies the command (removing cheetah comments, standardising galaxy dynamic vars)
    resolves aliases (temporary variables) back to original params
    """
    xmltool = load_xmltool(runtime.tool.tool_path)
    text = xmltool.raw_command
    text = simplify_cmd(text, 'main_statement')
    text = mark_main_statement(text, xmltool)
    text = simplify_cmd(text, 'parsing')
    text = resolve_aliases(text)
    return text

def load_templated_command_str(inputs_dict: dict[str, Any]) -> str:
    """
    loads <command> section of the active tool's XML for analysis.
    as above, except performs cheetah eval to simplify command.
    omits some simplification steps (e.g. cheetah comments) as many of these are handled during templating. 
    """
    xmltool = load_xmltool(runtime.tool.tool_path)
    text = xmltool.raw_command
    text = simplify_cmd(text, 'main_statement')
    text = mark_main_statement(text, xmltool)
    text = simplify_cmd(text, 'templating')
    text = sectional_evaluate(text, inputs=inputs_dict)
    text = simplify_cmd(text, 'parsing')
    text = resolve_aliases(text)
    return text

