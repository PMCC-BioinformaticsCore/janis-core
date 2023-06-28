

from janis_core.ingestion.galaxy import expressions
from janis_core.ingestion.galaxy.expressions.patterns import COMPOUND_OPT

from ..tokens import Token
from ..tokens import TokenType
from ...gxtool.param.InputParam import BoolParam, SelectParam


NON_VALUE_TOKENTYPES = set([
    TokenType.FUNCTION_CALL, 
    TokenType.BACKTICK_SHELL_STATEMENT, 
    TokenType.LINUX_TEE, 
    TokenType.LINUX_REDIRECT,
    TokenType.LINUX_STREAM_MERGE,
    TokenType.END_STATEMENT,
    TokenType.EXCISION,
])

# general utils

def is_bool_select(token: Token) -> bool:
    if token.ttype == TokenType.GX_INPUT:
        match token.gxparam:
            case BoolParam():
                return True
            case SelectParam():
                if len(token.gxparam.options) > 0:
                    return True
            case _:
                pass
    return False

def different_params(ctoken: Token, ntoken: Token) -> bool:
    # for when boolparam followed by a different gxparam
    if ctoken.gxparam and ntoken.gxparam:   # curr and next both have gxparams
        if type(ctoken.gxparam) != type(ntoken.gxparam):  # not same gxparam
            return True
    return False


# FLAG STUFF

def boolparam_then_gxparam(ctoken: Token, ntoken: Token) -> bool:
    # for when boolparam followed by a different gxparam
    if looks_like_a_flag(ctoken):
        if different_params(ctoken, ntoken):
            if isinstance(ctoken.gxparam, BoolParam):
                return True
    return False

def flag_then_flag(ctoken: Token, ntoken: Token) -> bool:
    if looks_like_a_flag(ctoken) and looks_like_a_flag(ntoken):
        return True
    return False

def flag_then_null(ctoken: Token, ntoken: Token) -> bool:
    if looks_like_a_flag(ctoken) and ntoken.ttype in NON_VALUE_TOKENTYPES:
        return True
    return False

def looks_like_a_flag(token: Token) -> bool:
    allowed_prefix_types = [TokenType.STRING, TokenType.INTEGER]
    if token.ttype == TokenType.FORCED_PREFIX:
        return True
    if token.ttype in allowed_prefix_types and token.text.startswith('-'):
        return True
    return False

flag_conditions = [
    boolparam_then_gxparam,
    flag_then_flag,
    flag_then_null,
]

def is_flag(ctoken: Token, ntoken: Token) -> bool:
    for condition in flag_conditions:
        if condition(ctoken, ntoken):
            return True
    return False


# OPTION STUFF 

def kvlinker(ctoken: Token, ntoken: Token) -> bool:
    if ntoken.ttype == TokenType.KV_LINKER:
        return True 
    return False

def compound_option(ctoken: Token, ntoken: Token) -> bool:
    if looks_like_a_flag(ctoken):
        if has_compound_structure(ctoken) and not is_positional(ntoken):
            return True
    return False

def has_compound_structure(token: Token) -> bool:
    compound_opts = expressions.get_matches(token.text, COMPOUND_OPT)
    if compound_opts:
        match = compound_opts[0]
        value = match.group(2)
        if int(value) > 3:
            return True
    return False

def prefix_then_value(ctoken: Token, ntoken: Token) -> bool:
    if not is_flag(ctoken, ntoken):
        if looks_like_a_flag(ctoken) and is_positional(ntoken):
            return True
    return False

option_conditions = [
    kvlinker,
    compound_option,
    prefix_then_value
]

def is_option(ctoken: Token, ntoken: Token) -> bool:
    """
    idetifies whether the current and next token implies an option.
    want to make sure it follows --prefix value structure (or variations).
    want to ensure other conditions arent True.
        - must be in same text construct (normal text, conditional block, loop)
        - must not be different galaxy params! this implies different tool args. 
    """
    if within_different_constructs(ctoken, ntoken):
        return False
    if different_params(ctoken, ntoken):
        return False
    for condition in option_conditions:
        if condition(ctoken, ntoken):
            return True
    return False

def within_different_constructs(ctoken: Token, ntoken: Token) -> bool:
    # one token is in a construct, other is not
    if ctoken.construct and not ntoken.construct:
        return True
    elif not ctoken.construct and ntoken.construct:
        return True
    # both tokens are in constructs, but the construct is not the same
    elif ctoken.construct and ntoken.construct:
        if ctoken.construct.uuid != ntoken.construct.uuid:
            return True
    return False


# POSITIONAL STUFF

def is_positional(token: Token) -> bool:
    if not looks_like_a_flag(token):
        if token.ttype not in NON_VALUE_TOKENTYPES:
            return True
    return False

