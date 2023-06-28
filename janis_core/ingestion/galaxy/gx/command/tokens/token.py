

import regex as re

from dataclasses import dataclass
from enum import Enum, auto
from typing import Optional

from janis_core.ingestion.galaxy.gx.gxtool.param import Param
from janis_core.ingestion.galaxy.gx.command.cmdstr import constructs

class TokenType(Enum):
    FUNCTION_CALL   = auto()
    SCRIPT          = auto()
    GX_INPUT        = auto()
    GX_OUTPUT       = auto()
    GX_KW_DYNAMIC   = auto()
    GX_KW_STATIC    = auto()
    ENV_VAR         = auto()
    KV_LINKER       = auto()
    FORCED_PREFIX   = auto()
    BACKTICK_SHELL_STATEMENT = auto()
    STRING          = auto()
    INTEGER         = auto()
    FLOAT           = auto()
    LINUX_TEE       = auto()
    LINUX_REDIRECT  = auto()
    LINUX_STREAM_MERGE = auto()
    START_STATEMENT = auto()
    EXCISION        = auto()
    END_STATEMENT   = auto()
    EMPTY_STRING    = auto()
    UNKNOWN         = auto()


@dataclass
class Token:
    match: re.Match[str]
    ttype: TokenType

    quoted: bool=False
    gxparam: Optional[Param]=None
    position: Optional[int] = None
    in_conditional: bool = False
    in_loop: bool = False
    construct: Optional[constructs.Construct] = None

    @property
    def text(self) -> str:
        return self.match[0] # type: ignore
    
    @property
    def start(self) -> int:
        return self.match.start()

    @property
    def end(self) -> int:
        return self.match.end()

    def __repr__(self) -> str:
        return f'Token: {self.text}, gxparam: {self.gxparam}'

