
from typing import Dict, List, Any, Optional, Union

from janis_core.operators.logical import If
from janis_core.types import Filename

from janis_core import (
    ToolOutput,
    ToolInput,
    CommandTool,
    InputSelector,
    WildcardSelector,
    StringFormatter,
    ToolArgument,
    InputDocumentation,
    InputQualityType,
)

from janis_core.types import (
    Stdout
)

class MessagingTestTool(CommandTool):
    def tool(self) -> str:
        return "MessagingTestTool"

    def base_command(self) -> Optional[Union[str, List[str]]]:
        return "echo"

    def arguments(self) -> List[ToolArgument]:
        return [ToolArgument("arg", position=0)]
    
    def inputs(self) -> List[ToolInput]:
        return [ToolInput("inp", str, position=0)]

    def outputs(self):
        return [ToolOutput("out", Stdout)]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"
    
mytool = MessagingTestTool()

print()