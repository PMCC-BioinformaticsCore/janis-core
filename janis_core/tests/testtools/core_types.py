
from typing import List, Optional, Union

from janis_core import (
    ToolOutput,
    ToolInput,
    CommandTool,
)

from janis_core.types import (
    File,
    String,
    Int,
    Stdout,
)


class TypeTestTool(CommandTool):
    def tool(self) -> str:
        return "TypeTestTool"

    def base_command(self) -> Optional[Union[str, List[str]]]:
        return "echo"

    def inputs(self) -> List[ToolInput]:
        return []

    def outputs(self):
        return [ToolOutput("out", Stdout)]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"


class FileTypeTestTool(TypeTestTool):
    def tool(self) -> str:
        return "FileTypeTestTool"

    def inputs(self) -> List[ToolInput]:
        return [ToolInput("inp", File, position=1)]


class StringTypeTestTool(TypeTestTool):
    def tool(self) -> str:
        return "StringTypeTestTool"

    def inputs(self) -> List[ToolInput]:
        return [ToolInput("inp", String, position=1)]


class IntTypeTestTool(TypeTestTool):
    def tool(self) -> str:
        return "IntTypeTestTool"

    def inputs(self) -> List[ToolInput]:
        return [ToolInput("inp", Int, position=1)]
