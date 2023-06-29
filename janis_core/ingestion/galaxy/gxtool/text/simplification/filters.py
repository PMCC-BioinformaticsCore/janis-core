

from .. import utils
import regex as re

from janis_core.ingestion.galaxy import expressions
from janis_core.ingestion.galaxy.expressions.patterns import (
    VARIABLES_FMT2,
    FUNCTION_CALL_FMT1,
    FUNCTION_CALL_FMT2,
    BACKTICK_SECTION,
    QUOTED_SECTION,
    GX_DYNAMIC_KEYWORDS,
    # GX_STATIC_KEYWORDS,
)

def flatten_nesting(cmdstr: str) -> str:
    cmdlines = cmdstr.split('\n')
    cmdlines = [ln.strip() for ln in cmdlines]
    cmdlines = [ln for ln in cmdlines if ln != '']
    return utils.join_lines(cmdlines)

def replace_function_calls(cmdstr: str) -> str:
    cmdlines = utils.split_lines(cmdstr)
    out: list[str] = []
    for line in cmdlines:
        matches = expressions.get_matches(line, FUNCTION_CALL_FMT1)
        matches += expressions.get_matches(line, FUNCTION_CALL_FMT2)
        for match in matches:
            # logging.has_cheetah_function()
            old_section = match[0]
            new_section = '__FUNCTION_CALL__'
            line = line.replace(old_section, new_section)
        out.append(line)
    return utils.join_lines(out)

def replace_backticks(cmdstr: str) -> str:
    matches = expressions.get_matches(cmdstr, BACKTICK_SECTION)
    for match in matches:
        # logging.has_backtick_statement()
        old_section = match[0]
        new_section = '__BACKTICK_SHELL_STATEMENT__'
        cmdstr = cmdstr.replace(old_section, new_section)
    return cmdstr

def interpret_raw(cmdstr: str) -> str:
    return cmdstr.replace('\\', '')

def flatten_multiline_strings(cmdstr: str) -> str:
    matches = expressions.get_matches(cmdstr, QUOTED_SECTION)
    for match in matches:
        if '\n' in match[0]:
            # logging.has_multiline_str()
            old_section = match[0]
            new_section = match[0].replace('\n', ' ')
            cmdstr = cmdstr.replace(old_section, new_section)
    return cmdstr

def translate_variable_markers(cmdstr: str) -> str:
    return cmdstr.replace("gxparam_", "$")

def standardise_variable_format(cmdstr: str) -> str:
    """
    modifies cmd word to ensure the $var format is present, rather than ${var}
    takes a safe approach using regex and resolving all vars one by one
    """
    iterator = re.finditer(VARIABLES_FMT2, cmdstr)
    matches = [m for m in iterator]
    matches = sorted(matches, key=lambda m: m.start(), reverse=True)

    for match in matches:
        new_segment = f'${match.group(1)}'
        cmdstr = cmdstr[:match.start()] + new_segment + cmdstr[match.end():]
    return cmdstr

def remove_empty_quotes(cmdstr: str) -> str:
    cmdstr = cmdstr.replace('""', '')
    cmdstr = cmdstr.replace("''", '')
    return cmdstr

def simplify_sh_constructs(cmdstr: str) -> str:
    """
    this function standardises the different equivalent 
    forms of linux operations into a single common form
    """
    cmdstr = cmdstr.replace("&amp;", "&")
    cmdstr = cmdstr.replace("&lt;", "<")
    cmdstr = cmdstr.replace("&gt;", ">")
    cmdstr = cmdstr.replace("|&", "2>&1 |")
    cmdstr = cmdstr.replace("| tee", "|tee")
    cmdstr = cmdstr.replace("1>", ">")
    return cmdstr 

def simplify_galaxy_dynamic_vars(cmdstr: str) -> str:
    """  ${GALAXY_SLOTS:-2} -> 2   etc """
    matches = expressions.get_matches(cmdstr, GX_DYNAMIC_KEYWORDS)
    for match in matches:
        cmdstr = cmdstr.replace(match[0], match.group(1)) 
    return cmdstr

def remove_cheetah_comments(cmdstr: str) -> str:
    """
    removes cheetah comments from shellparser lines
    comments can be whole line, or part way through
    """
    cmdlines: list[str] = utils.split_lines(cmdstr)
    outlines: list[str] = []

    for line in cmdlines:
        comment_start, _ = expressions.find_unquoted(line, '##')
        if comment_start != -1:
            # override line with comment removed
            line = line[:comment_start].strip()
        # make sure we didnt trim a full line comment and now its an empty string
        if line != '':
            outlines.append(line)
    return utils.join_lines(outlines)


