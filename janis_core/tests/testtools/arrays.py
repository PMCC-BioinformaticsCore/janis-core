


from typing import Optional 

from janis_core import (
    ToolOutput,
    ToolInput,
    CommandTool,
    WildcardSelector
)

from janis_core.types import (
    File,
    String,
    Int,
    Stdout,
    Array
)


class InArrayTestTool(CommandTool):
    def tool(self) -> str:
        return "InArrayTestTool"

    def base_command(self) -> Optional[str | list[str]]:
        return "echo"

    def inputs(self) -> list[ToolInput]:
        return [ToolInput("in", Array(File))]

    def outputs(self):
        return [ToolOutput("out", Stdout)]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"


class OutArrayTestTool(CommandTool):
    def tool(self) -> str:
        return "OutArrayTestTool"

    def base_command(self) -> Optional[str | list[str]]:
        return "echo"

    def inputs(self) -> list[ToolInput]:
        return [ToolInput("in", Array(File))]

    def outputs(self):
        return [ToolOutput("out", Array(String), selector=WildcardSelector('*'))]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"


