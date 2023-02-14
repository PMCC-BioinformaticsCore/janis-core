

from janis_core.ingestion.galaxy.gx.gxtool import XMLToolDefinition
from janis_core.ingestion.galaxy.gx.gxtool.param import ParamRegister
from janis_core.ingestion.galaxy.gx.gxtool.TestRegister import TestRegister
from .mock_metadata import MOCK_TOOL_METADATA


MOCK_XMLTOOL = XMLToolDefinition(
    metadata=MOCK_TOOL_METADATA,
    raw_command='',
    configfiles=[],
    inputs=ParamRegister(),
    outputs=ParamRegister(),
    tests=TestRegister([])
)