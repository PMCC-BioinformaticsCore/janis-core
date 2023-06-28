


from .token import Token, TokenType

from janis_core.ingestion.galaxy import expressions
from janis_core.ingestion.galaxy.expressions.patterns import ALL


def spawn_end_sentinel() -> Token:
    matches = expressions.get_matches('end', ALL)
    return Token(matches[0], TokenType.END_STATEMENT)

def spawn_kv_linker(delim: str) -> Token:
    matches = expressions.get_matches(delim, ALL)
    return Token(matches[0], TokenType.KV_LINKER)

def spawn_function_call() -> Token:
    matches = expressions.get_matches('__FUNCTION_CALL__', ALL)
    return Token(matches[0], TokenType.FUNCTION_CALL)
    
def spawn_backtick_section() -> Token:
    matches = expressions.get_matches('__BACKTICK_SHELL_STATEMENT__', ALL)
    return Token(matches[0], TokenType.BACKTICK_SHELL_STATEMENT)

