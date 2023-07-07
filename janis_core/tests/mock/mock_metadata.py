

from janis_core.ingestion.galaxy.gxtool.model import XMLMetadata
from janis_core.ingestion.galaxy.internal_model.workflow import WorkflowMetadata


MOCK_WORKFLOW_METADATA = WorkflowMetadata(
    name='assembly',
    uuid='1',
    version='1.1',
    annotation='a hybrid assembly workflow using unicycler',
    tags=['assembly']
)

MOCK_TOOL_METADATA = XMLMetadata(
    id='abricate',
    name='ABRicate',
    version='1.0.1',
    description='Mass screening of contigs for antimicrobial and virulence genes',
    help='\n**abricate help**\n ',
    requirements=[],
    citations=[],
    creator=None
)