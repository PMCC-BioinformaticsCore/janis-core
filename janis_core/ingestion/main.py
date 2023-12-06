
from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from janis_core import Tool

def ingest_galaxy(uri: str) -> Tool:
    from .galaxy import parse_galaxy
    return parse_galaxy(uri)

def ingest_cwl(path: str) -> Tool:
    from .cwl import parse as parse_cwl
    return parse_cwl(path)

def ingest_wdl(path: str) -> Tool:
    from .wdl import WdlParser
    return WdlParser.from_doc(path)

def ingest(
    path: str, 
    fmt: str, 
    build_galaxy_tool_images: bool = False, 
    ) -> Tool:
    from janis_core import settings
    from janis_core.messages import configure_logging
    from .SupportedIngestion import SupportedIngestion


    # setup logging
    configure_logging()                         
    
    # set ingest settings
    settings.ingest.SOURCE = fmt                     
    settings.validation.STRICT_IDENTIFIERS = False
    settings.validation.VALIDATE_STRINGFORMATTERS = False
    if build_galaxy_tool_images:
        settings.ingest.galaxy.GEN_IMAGES = True

    # do ingest
    assert(fmt in SupportedIngestion.all())  # validate format

    if fmt == 'galaxy':
        return ingest_galaxy(path)
    elif fmt == 'cwl':
        return ingest_cwl(path)
    elif fmt == 'wdl':
        return ingest_wdl(path)
    raise Exception