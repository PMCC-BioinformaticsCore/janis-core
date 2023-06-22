

from typing import Optional
from janis_core.ingestion.galaxy import logging
from janis_core.ingestion.galaxy import runtime
from janis_core.ingestion.galaxy.utils import galaxy as utils


def tool_setup(args: dict[str, Optional[str]]) -> None:
    runtime.tool.set(from_args=args)
    logging.msg_parsing_tool()
    validate_tool_settings()


### VALIDATION ###

def validate_tool_settings() -> None:
    # both local & remote params not given
    if not _has_xml() or not _valid_xml():
        raise RuntimeError('no valid xml file')

def _has_xml() -> bool:
    if runtime.tool.tool_path:
        return True
    return False

def _valid_xml() -> bool:
    path = runtime.tool.tool_path
    if utils.is_tool_xml(path):
        return True
    return False


