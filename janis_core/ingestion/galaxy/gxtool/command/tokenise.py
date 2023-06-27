

from typing import Optional

from janis_core.ingestion.galaxy.gxtool.model import XMLTool
from .cmdstr.RealisedTokenValues import RealisedTokens
from .cmdstr.RealisedTokenValues import RealisedTokenFactory


def tokenise_text(text: str, xmltool: Optional[XMLTool]=None) -> list[RealisedTokens]:
    factory = RealisedTokenFactory(xmltool)
    return factory.try_tokenise_text(text)

def tokenise_line(line: str, xmltool: Optional[XMLTool]=None) -> list[RealisedTokens]:
    factory = RealisedTokenFactory(xmltool)
    return factory.try_tokenise_line(line)