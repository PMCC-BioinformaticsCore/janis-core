

from galaxy2janis.gx.gxtool import XMLToolDefinition
from galaxy2janis.gx.gxtool.param import ParamRegister
from galaxy2janis.gx.gxtool.TestRegister import TestRegister
from .mock_metadata import MOCK_TOOL_METADATA


MOCK_XMLTOOL = XMLToolDefinition(
    metadata=MOCK_TOOL_METADATA,
    raw_command='',
    configfiles=[],
    inputs=ParamRegister(),
    outputs=ParamRegister(),
    tests=TestRegister([])
)