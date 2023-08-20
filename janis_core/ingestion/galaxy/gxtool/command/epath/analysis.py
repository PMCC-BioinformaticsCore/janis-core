
from typing import Tuple, Optional
from copy import deepcopy
from janis_core.ingestion.galaxy import expressions

from ...model import XMLBoolParam, XMLSelectParam
from ..tokens import Token
from ..tokens import TokenType
from ..Command import Command
from .ExecutionPath import EPathPosition


NON_VALUE_TOKENTYPES = set([
    TokenType.END_STATEMENT,
    TokenType.LINUX_TEE, 
    TokenType.LINUX_REDIRECT,
    TokenType.LINUX_STREAM_MERGE,
    TokenType.EXCISION,
])

STATEMENT_END_TOKENTYPES = set([
    TokenType.END_STATEMENT,
    TokenType.LINUX_TEE, 
    TokenType.LINUX_REDIRECT,
    TokenType.LINUX_STREAM_MERGE,
])

"""
utility functions to analyse an ExecutionPath. 
this exists so we can extract components from a list of tokens. 

includes functionality to:
- detect possible component types
- identify which tokens relate to previously extracted components

"""


### --- GENERAL UTILS --- ###

# attributes

def different_params(cpos: EPathPosition, npos: EPathPosition) -> bool:
    """detect when current & next have gxparams, but they're different"""
    if cpos.token.gxparam and npos.token.gxparam:   # curr and next both have gxparams
        if type(cpos.token.gxparam) != type(npos.token.gxparam):  # not same gxparam
            return True
    return False

def different_components(cpos: EPathPosition, npos: EPathPosition) -> bool:
    """detect when current & next have different components"""
    # this one no component, next has component
    if not cpos.component and npos.component:
        return True
    # this one has component, next has no component
    elif cpos.component and not npos.component:
        return True
    # both have components, but theyre different
    elif cpos.component and npos.component and cpos.component.uuid != npos.component.uuid:
        return True
    return False

def different_constructs(cpos: EPathPosition, npos: EPathPosition) -> bool:
    """detect when current & next exist in different cheetah constructs"""
    # one token is in a construct, other is not
    if cpos.token.construct and not npos.token.construct:
        return True
    elif not cpos.token.construct and npos.token.construct:
        return True
    # both pos.tokens are in constructs, but the construct is not the same
    elif cpos.token.construct and npos.token.construct:
        if cpos.token.construct.uuid != npos.token.construct.uuid:
            return True
    return False

def different_lines(cpos: EPathPosition, npos: EPathPosition) -> bool:
    """detect when current & next exist in different cheetah constructs"""
    if cpos.token.line != npos.token.line:
        return True
    return False

# token structure

def is_prefix(pos: EPathPosition) -> bool:
    allowed_prefix_types = [TokenType.STRING, TokenType.INTEGER]
    if pos.token.ttype == TokenType.FORCED_PREFIX:
        return True
    if pos.token.ttype in allowed_prefix_types and pos.token.text.startswith('-'):
        return True
    return False

def has_kvlinker_structure(cpos: EPathPosition, npos: EPathPosition) -> bool:
    if npos.token.ttype == TokenType.KV_LINKER:
        if cpos.token.ttype == TokenType.FORCED_PREFIX:
            return True
        elif cpos.token.ttype == TokenType.STRING:
            return True
        elif cpos.token.ttype == TokenType.INTEGER and cpos.token.text.startswith('-'):
            return True
    return False

def has_prefix_value_structure(cpos: EPathPosition, npos: EPathPosition) -> bool:
    if is_prefix(cpos) and is_positional(npos):
        return True
    return False

def has_compound_structure(cpos: EPathPosition, npos: EPathPosition) -> bool:
    if is_prefix(cpos) and not is_positional(npos):
        compound_opts = expressions.get_matches(cpos.token.text, expressions.patterns.COMPOUND_OPT)
        if compound_opts:
            match = compound_opts[0]
            value = match.group(2)
            if int(value) > 3:
                return True
    return False


### --- POSITIONALS --- ###

# identifying new positional components

def is_positional(pos: EPathPosition) -> bool:
    if is_prefix(pos):
        return False
    if pos.token.ttype in NON_VALUE_TOKENTYPES:
        return False
    return True


### --- FLAGS --- ###

# identifying new flag components

def is_flag(cpos: EPathPosition, npos: EPathPosition) -> bool:
    """identifies whether the current and next position implies a flag."""
    # current must be a prefix
    if not is_prefix(cpos):
        return False
    
    # positive indicators
    if is_prefix(npos):
        return True
    elif npos.token.ttype in NON_VALUE_TOKENTYPES:
        return True
    elif different_components(cpos, npos):
        return True
    elif different_constructs(cpos, npos):
        return True
    elif different_params(cpos, npos):
        return True
    
    return False


### --- OPTIONS --- ###

# identifying new option components

def is_option(cpos: EPathPosition, npos: EPathPosition) -> bool:
    """
    idetifies whether the current and next token implies an option.
    want to make sure it follows --prefix value structure (or variations).
    want to ensure other conditions arent True.
        - must be in same text construct (normal text, conditional block, loop)
        - must not be different galaxy params! this implies different tool args. 
    """
    # negative indicators
    if different_constructs(cpos, npos):
        return False
    if different_components(cpos, npos):
        return False
    if different_params(cpos, npos):
        return False
    
    # positive indicators
    if has_prefix_value_structure(cpos, npos):
        return True
    if has_kvlinker_structure(cpos, npos):
        return True
    if has_compound_structure(cpos, npos):
        return True
    
    return False

def get_option_values_span(ptr: int, positions: list[EPathPosition]) -> Tuple[int, int]:
    cpos = positions[ptr]
    npos = positions[ptr + 1]

    if has_kvlinker_structure(cpos, npos):
        start, stop = ptr + 2, ptr + 2
    elif has_prefix_value_structure(cpos, npos):
        start, stop = get_option_values_span_greedy(ptr + 1, positions)
    elif has_compound_structure(cpos, npos):
        start, stop = ptr, ptr
    else:
        raise RuntimeError('should not be here')
    
    return start, stop

def get_option_values_span_greedy(ptr: int, positions: list[EPathPosition]) -> Tuple[int, int]:
    start = deepcopy(ptr)
    while should_eat_value(ptr, ptr + 1, positions):
        ptr += 1
    return start, ptr

def should_eat_value(curr_ptr: int, next_ptr: int, positions: list[EPathPosition]) -> bool:
    if not is_positional(positions[next_ptr]):
        return False
    if different_constructs(positions[curr_ptr], positions[next_ptr]):
        return False
    if different_components(positions[curr_ptr], positions[next_ptr]):
        return False
    if different_lines(positions[curr_ptr], positions[next_ptr]):
        return False
    if is_last_positional(next_ptr, positions):
        return False
    return True

def is_last_positional(ptr: int, positions: list[EPathPosition]) -> bool:
    if is_positional(positions[ptr]):
        npos = positions[ptr + 1]
        if npos.token.ttype in STATEMENT_END_TOKENTYPES:
            return True
    return False




