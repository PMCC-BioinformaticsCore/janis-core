

from galaxy2janis.gx.gxtool import ToolXMLMetadata
from galaxy2janis.entities.workflow import WorkflowMetadata


MOCK_WORKFLOW_METADATA = WorkflowMetadata(
    name='assembly',
    uuid='1',
    version='1.1',
    annotation='a hybrid assembly workflow using unicycler',
    tags=['assembly']
)

MOCK_TOOL_METADATA = ToolXMLMetadata(
    id='abricate',
    name='ABRicate',
    version='1.0.1',
    description='Mass screening of contigs for antimicrobial and virulence genes',
    help='\n**abricate help**\n ',
    requirements=[],
    citations=[],
    creator=None
)