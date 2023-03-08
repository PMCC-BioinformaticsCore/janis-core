
import os 
from typing import Optional

from janis_core import Tool
from janis_core.messages import configure_logging
from janis_core.messages import info_ingesting_tool
from janis_core.messages import info_ingesting_workflow

from janis_core.ingestion.galaxy import ingest_tool
from janis_core.ingestion.galaxy import ingest_workflow
from janis_core.ingestion.galaxy.janis_mapping import to_janis_tool
from janis_core.ingestion.galaxy.janis_mapping import to_janis_workflow
from janis_core import settings

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
    # strict_identifiers: Optional[bool]=False,
    # allow_incorrect_number_of_sources: Optional[bool]=True,
    # allow_non_array_scatter_input: Optional[bool]=True
) -> Tool:
    
    # settings
    settings.identifiers.STRICT_IDENTIFIERS = False
    settings.graph.ALLOW_INCORRECT_NUMBER_OF_SOURCES = True
    settings.graph.ALLOW_NON_ARRAY_SCATTER_INPUT = True
    settings.graph.ALLOW_INCOMPATIBLE_TYPES = True

    DEV_MODE = True
    if DEV_MODE:
        settings.ingest.cwl.INGEST_JAVASCRIPT_EXPRESSIONS = True
        settings.ingest.cwl.REQUIRE_CWL_VERSION = False
        settings.datatypes.ALLOW_UNPARSEABLE_DATATYPES = True
        settings.graph.ALLOW_UNKNOWN_SOURCE = True
        settings.graph.ALLOW_UNKNOWN_SCATTER_FIELDS = True
    else:
        settings.ingest.cwl.INGEST_JAVASCRIPT_EXPRESSIONS = False
        settings.ingest.cwl.REQUIRE_CWL_VERSION = False
        settings.datatypes.ALLOW_UNPARSEABLE_DATATYPES = True
        settings.graph.ALLOW_UNKNOWN_SOURCE = False
        settings.graph.ALLOW_UNKNOWN_SCATTER_FIELDS = False
        
    # setup
    configure_logging()
    
    # validation
    assert(format in SupportedIngestion.all())  # validate format

    # ingest
    ingest_func = ingestor_map[format]          # select ingestor
    internal = ingest_func(path)                # ingest
    return internal