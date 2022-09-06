


from galaxy2janis.entities.tool import Tool

from .mock_params import MOCK_PARAM_REGISTER
from .mock_metadata import MOCK_TOOL_METADATA
from .mock_components import (
    MOCK_POSITIONAL1,
    MOCK_FLAG1,
    MOCK_OPTION1,
    MOCK_OPTION2,
    MOCK_REDIRECT_OUTPUT,
    MOCK_WILDCARD_OUTPUT,
    MOCK_INPUT_OUTPUT
)


MOCK_TOOL_ABRICATE = Tool(
    metadata=MOCK_TOOL_METADATA,
    gxparam_register=MOCK_PARAM_REGISTER,
    configfiles=[],
    container='quay.io/biocontainers/abricate:1.0.1--ha8f3691_1',
    base_command=['abricate'],
)
# dont change this order
MOCK_TOOL_ABRICATE.add_input(MOCK_POSITIONAL1)
MOCK_TOOL_ABRICATE.add_input(MOCK_FLAG1)
MOCK_TOOL_ABRICATE.add_input(MOCK_OPTION1)
MOCK_TOOL_ABRICATE.add_input(MOCK_OPTION2)
MOCK_TOOL_ABRICATE.add_output(MOCK_REDIRECT_OUTPUT)
MOCK_TOOL_ABRICATE.add_output(MOCK_WILDCARD_OUTPUT)
MOCK_TOOL_ABRICATE.add_output(MOCK_INPUT_OUTPUT)