
from janis_core import CommandToolBuilder, ToolInput, String

MOCK_TOOL1 = CommandToolBuilder(
    tool='mock',
    base_command=["mock"],
    inputs=[ToolInput("inp1", String())], 
    outputs=[],
    container='ubuntu:latest',
    version='0.1'
)

