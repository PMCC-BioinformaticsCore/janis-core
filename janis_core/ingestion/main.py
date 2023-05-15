
import os 

from janis_core import Tool
from janis_core.messages import configure_logging
from janis_core.messages import info_ingesting_tool
from janis_core.messages import info_ingesting_workflow

from janis_core.ingestion.galaxy import ingest_tool
from janis_core.ingestion.galaxy import ingest_workflow
from janis_core.ingestion.galaxy.janis_mapping import to_janis_tool
from janis_core.ingestion.galaxy.janis_mapping import to_janis_workflow
 
from .SupportedIngestion import SupportedIngestion
from .cwl import parse as parse_cwl
from .wdl import WdlParser

def ingest_galaxy(path: str) -> Tool:
    # this function is essentially the same as CWlParser / WdlParser
    name, ext = os.path.splitext(path)  # is this a tool .xml, or a .ga workflow? 
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

def ingest_cwl(path: str) -> Tool:
    return parse_cwl(path)

def ingest_wdl(path: str) -> Tool:
    return WdlParser.from_doc(path)


# this is using the SupportedIngestion Enum values. 
# should it just be the Enum? these strings are in 2 places now.  
ingestor_map = {  
    'galaxy': ingest_galaxy,
    'cwl': ingest_cwl,
    'wdl': ingest_wdl
}


def ingest(
    path: str, 
    format: str,
) -> Tool:
    
    # setup
    configure_logging()
    
    # validation
    assert(format in SupportedIngestion.all())  # validate format

    # ingest
    ingest_func = ingestor_map[format]          # select ingestor
    internal = ingest_func(path)                # ingest
    return internal