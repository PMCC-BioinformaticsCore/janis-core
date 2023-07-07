

from janis_core.ingestion.galaxy.gxtool.model import XMLTool
from janis_core.ingestion.galaxy.gxtool.model import XMLParamRegister
from janis_core.ingestion.galaxy.gxtool.model import XMLTestRegister
from .mock_metadata import MOCK_TOOL_METADATA


MOCK_XMLTOOL = XMLTool(
    metadata=MOCK_TOOL_METADATA,
    raw_command='',
    configfiles=[],
    scripts=[],
    inputs=XMLParamRegister(),
    outputs=XMLParamRegister(),
    tests=XMLTestRegister([])
)