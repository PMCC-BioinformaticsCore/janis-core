
import os 

from janis_core import CommandTool
from janis_core import Workflow
from janis_core.messages import configure_logging
from janis_core.messages import info_ingesting_tool
from janis_core.messages import info_ingesting_workflow

from galaxy2janis import ingest_tool
from galaxy2janis import ingest_workflow

from .SupportedIngestion import SupportedIngestion
from .janis_mapping import to_janis_tool
from .janis_mapping import to_janis_workflow
from .fromcwl import CWlParser
from .fromwdl import WdlParser


def ingest_galaxy(path: str) -> CommandTool | Workflow:
    # this function is essentially the same as CWlParser / WdlParser
    name, ext = os.path.splitext(path)
    if ext == '.xml':
        info_ingesting_tool('galaxy', name)
        tool = ingest_tool(path)
        return to_janis_tool(tool)
    elif ext == '.ga':
        info_ingesting_workflow('galaxy', name)
        workflow = ingest_workflow(path)
        return to_janis_workflow(workflow)
    else:
        raise ValueError("file must end in '.xml' or '.ga' for galaxy ingestion")

def ingest_cwl(path: str) -> CommandTool | Workflow:
    return CWlParser.from_doc(path)

def ingest_wdl(path: str) -> CommandTool | Workflow:
    return WdlParser.from_doc(path)


# this is using the SupportedIngestion Enum values. 
# should it just be the Enum? these strings are in 2 places now.  
ingestor_map = {  
    'galaxy': ingest_galaxy,
    'cwl': ingest_cwl,
    'wdl': ingest_wdl
}


def ingest(path: str, format: str) -> CommandTool | Workflow:
    configure_logging()
    assert(format in SupportedIngestion.all())  # validate format
    ingest_func = ingestor_map[format]          # select ingestor
    internal = ingest_func(path)                # ingest
    return internal