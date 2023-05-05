

from typing import Optional, Any
from janis_core.ingestion.galaxy.gx.gxtool import XMLToolDefinition
from janis_core.ingestion.galaxy.gx.command import Command

# this module imports
from .Tool import Tool
from .ToolFactory import ToolFactory


def gen_tool(
    xmltool: XMLToolDefinition, 
    command: Command, 
    container: Optional[str], 
    gxstep: Optional[dict[str, Any]]=None
    ) -> Tool:
    factory = ToolFactory(xmltool, command, container, gxstep)
    tool = factory.create()
    return tool

