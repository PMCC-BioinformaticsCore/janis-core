

from __future__ import annotations
from typing import TYPE_CHECKING, Optional, Tuple
if TYPE_CHECKING:
    from janis_core.ingestion.galaxy.gxtool.model import XMLTool

import regex as re

from janis_core.ingestion.galaxy import expressions
from janis_core.ingestion.galaxy.expressions.patterns import (
    INTEGER,
    FLOAT,
    STRING,
    LINUX_REDIRECT,
    LINUX_TEE,
    LINUX_STREAM_MERGE,
    SCRIPT,
    GX_DYNAMIC_KEYWORDS,
    GX_STATIC_KEYWORDS,
    EMPTY_STRING,
    VARIABLES_FMT1,
    ALL,
)

from .token import Token, TokenType
from .prioritisation import TokenOrderingStrategy
from .prioritisation import FirstTokenOrderingStrategy
from .prioritisation import LongestTokenOrderingStrategy
from .prioritisation import PriorityTokenOrderingStrategy
from . import utils as token_utils


def tokenize(word: str, prioritisation: str='priority', xmltool: Optional[XMLTool]=None) -> list[Token]:
    quoted, word = _handle_quoted(word)
    tokens: list[Token] = []
    tokens += _spawn_kv_tokens(word, prioritisation, xmltool)
    if not tokens:
        token = _spawn_single(word, prioritisation, xmltool)
        tokens += [token]
    _mark_quoted(quoted, tokens)
    return tokens

def tokenize_single(word: str, prioritisation: str='priority', xmltool: Optional[XMLTool]=None) -> Token:
    quoted, word = _handle_quoted(word)
    token = _spawn_single(word, prioritisation)
    _mark_quoted(quoted, [token])
    return token

def _handle_quoted(word: str) -> Tuple[bool, str]:
    quotes = ['"', "'"]
    if word[0] in quotes and word[-1] == word[0]:
        word = word[1:-1]
        quoted = True
    else:
        quoted = False
    return quoted, word

def _mark_quoted(quoted: bool, tokens: list[Token]) -> None:
    if quoted:
        for token in tokens:
            token.quoted = True

generic_patterns = [
    (INTEGER, TokenType.INTEGER),
    (FLOAT, TokenType.FLOAT),
    (STRING, TokenType.STRING),
    (LINUX_REDIRECT, TokenType.LINUX_REDIRECT),
    (LINUX_TEE, TokenType.LINUX_TEE),
    (LINUX_STREAM_MERGE, TokenType.LINUX_STREAM_MERGE),
    (SCRIPT, TokenType.SCRIPT),
    (GX_DYNAMIC_KEYWORDS, TokenType.GX_KW_DYNAMIC),
    (GX_STATIC_KEYWORDS, TokenType.GX_KW_STATIC),
    (EMPTY_STRING, TokenType.EMPTY_STRING),
    (ALL, TokenType.UNKNOWN)
]

ordering_strategies: dict[str, TokenOrderingStrategy] = {
    'first':  FirstTokenOrderingStrategy(),
    'priority':  PriorityTokenOrderingStrategy(),
    'longest':  LongestTokenOrderingStrategy()
}

def _spawn_kv_tokens(word: str, prioritisation: str='priority', xmltool: Optional[XMLTool]=None) -> list[Token]:
    left_token_allowed_types = [
        TokenType.STRING,
        TokenType.INTEGER,
    ]
    right_token_allowed_types = [
        TokenType.STRING,
        TokenType.INTEGER,
        TokenType.FLOAT,
        TokenType.SCRIPT,
        TokenType.GX_INPUT,
        TokenType.GX_OUTPUT,
        TokenType.GX_KW_DYNAMIC,
        TokenType.GX_KW_STATIC,
        TokenType.ENV_VAR
    ]
    delims = ['=', ':']
    for delim in delims:
        if delim in word:
            left, right = word.split(delim, 1)
            left_token = _spawn_single(left, prioritisation, xmltool)
            right_token = _spawn_single(right, prioritisation, xmltool)
            if left_token.ttype in left_token_allowed_types:
                if right_token.ttype in right_token_allowed_types:
                    return [left_token, token_utils.spawn_kv_linker(delim), right_token]
    return []

def _spawn_single(word: str, prioritisation: str='priority', xmltool: Optional[XMLTool]=None) -> Token:
    """
    extracts the best token from a word.
    where multiple token types are possible, selection can be made 
    """
    token_list = _get_all_tokens(word, xmltool)
    token_list = _perform_default_ordering(token_list)
    final_ordering = _perform_final_ordering(token_list, prioritisation)
    return final_ordering[0]

def _get_all_tokens(the_string: str, xmltool: Optional[XMLTool]=None) -> list[Token]:
    """gets all the possible token interpretations of the_string"""  
    tokens: list[Token] = []
    tokens += _get_dunder_token(the_string)
    tokens += _get_generic_tokens(the_string)
    tokens += _get_variable_tokens(the_string, xmltool)
    return tokens

def _get_dunder_token(the_string: str) -> list[Token]:
    if the_string == '__FUNCTION_CALL__':
        return [token_utils.spawn_function_call()]
    if the_string == '__BACKTICK_SHELL_STATEMENT__':
        return [token_utils.spawn_backtick_section()]
    return []

def _get_generic_tokens(the_string: str) -> list[Token]:
    """gets all tokens except galaxy/env variables"""
    tokens: list[Token] = []
    for pattern, ttype in generic_patterns:
        matches = expressions.get_matches(the_string, pattern)
        tokens += [Token(m, ttype) for m in matches]
    return tokens

def _get_variable_tokens(the_string: str, xmltool: Optional[XMLTool]=None) -> list[Token]:
    """gets tokens for galaxy/env variables"""
    tokens: list[Token] = []
    matches = expressions.get_matches(the_string, VARIABLES_FMT1)
    base_vars = [strip_to_base_variable(m) for m in matches]
    base_vars = [v for v in base_vars if v is not None]
    
    for m, varname in zip(matches, base_vars):
        if xmltool:
            if xmltool.inputs.get(varname):
                tokens.append(Token(m, TokenType.GX_INPUT, gxparam=xmltool.inputs.get(varname)))
            elif xmltool.outputs.get(varname):
                tokens.append(Token(m, TokenType.GX_OUTPUT, gxparam=xmltool.outputs.get(varname)))
            else:
                tokens.append(Token(m, TokenType.ENV_VAR))
        else:
            tokens.append(Token(m, TokenType.ENV_VAR))
    return tokens

def strip_to_base_variable(match: re.Match[str]) -> Optional[str]:
    """trims function calls, attributes from variable matches"""
    text: str = match[0]
    text = _strip_quotes(text)
    text = _strip_braces(text)
    text = _strip_method_calls(text, match)
    text = _strip_common_attributes(text)
    text = text.replace('$', '')
    if text != '':
        return text
    return None

def _strip_quotes(text: str) -> str:
    text = text.replace('"', '')
    text = text.replace("'", '')
    return text

def _strip_braces(text: str) -> str:
    if text[1] == '{':
        text = text[0] + text[2:]
    if text[-1] == '}':
        text = text[:-1]
    # text = text.replace('{', '')
    # text = text.replace('}', '')
    return text

def _strip_method_calls(text: str, match: re.Match[str]) -> str:
    """
    only want cheetah variable references.  
    sometimes we see python string methods attached to a galaxy param var or cheetah functions (which have similar syntax to other vars of course). Want to remove these.
    """
    if match.end() < len(match.string) and match.string[match.end()] == '(':
        # object method?
        if '.' in text:
            # strip back method call
            text = text.rsplit('.', 1)[0]
        else:
            # is cheetah func call.  
            text = ''
    return text

def _perform_default_ordering(token_list: list[Token]) -> list[Token]:
    #default orderings (low to high priority) first, longest, priority
    for strategy in ordering_strategies.values():
        token_list = strategy.order(token_list)
    return token_list
    
def _perform_final_ordering(token_list: list[Token],  prioritisation: str) -> list[Token]:
    # overriding final prioritisation
    return ordering_strategies[prioritisation].order(token_list)

def _strip_common_attributes(text: str) -> str:
    #return text
    gx_attributes = set([
        '.forward',
        '.reverse',
        '.ext',
        '.value',
        '.name',
        #'.files_path',
        '.element_identifier'
    ])
    # needs to be recursive so we can iterately peel back 
    # eg  in1.forward.ext
    # need to peel .ext then peel .forward.
    for att in gx_attributes:
        if text.endswith(att):
            # strip from the right - num of chars in the att
            text = text[:-len(att)]
            # recurse
            text = _strip_common_attributes(text)
    return text
