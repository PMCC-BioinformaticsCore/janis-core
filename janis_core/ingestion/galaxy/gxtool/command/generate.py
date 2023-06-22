

from typing import Optional, Any
from janis_core.ingestion.galaxy.gxtool.model import XMLTool
from .factory import CommandFactory
from .Command import Command

 
def gen_command(xmltool: XMLTool, gxstep: Optional[dict[str, Any]]=None) -> Command:
    factory = CommandFactory(xmltool, gxstep)
    return factory.create()



