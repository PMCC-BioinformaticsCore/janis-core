
from typing import Optional

from janis_bioinformatics.data_types.bam import BamBai

from janis_core import (
    ToolOutput,
    ToolInput,
    CommandTool,
    ToolArgument,
    WildcardSelector,
    InputSelector

)

from .types import (
    SecondaryTestType,
    AppendedSecondaryTestType,
    ReplacedSecondaryTestType,
    NonEscapedSecondaryTestType
)

from janis_core.types import (
    File,
    String,
    Int,
    Stdout,
    Boolean,
    Filename
)


class TypeTestTool(CommandTool):
    def tool(self) -> str:
        return "TypeTestTool"

    def base_command(self) -> Optional[str | list[str]]:
        return "echo"

    def inputs(self) -> list[ToolInput]:
        return []

    def outputs(self):
        return [ToolOutput("out", Stdout)]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"


class FileTestTool(TypeTestTool):
    def tool(self) -> str:
        return "FileTypeTestTool"

    def inputs(self) -> list[ToolInput]:
        return [ToolInput("inp", File, position=1)]


class StringTestTool(TypeTestTool):
    def tool(self) -> str:
        return "StringTypeTestTool"

    def inputs(self) -> list[ToolInput]:
        return [ToolInput("inp", String, position=1)]


class IntTestTool(TypeTestTool):
    def tool(self) -> str:
        return "IntTypeTestTool"

    def inputs(self) -> list[ToolInput]:
        return [ToolInput("inp", Int, position=1)]


class ComponentsTestTool(CommandTool):
    def tool(self) -> str:
        return "ComponentTypeTestTool"

    def base_command(self) -> Optional[str | list[str]]:
        return "echo"

    def inputs(self) -> list[ToolInput]:
        return [
            ToolInput("pos_basic", File, position=1),
            ToolInput("pos_default", String, default="DEFAULT", position=2),
            ToolInput("pos_optional", String(optional=True), position=3),

            ToolInput("flag_true", Boolean, position=4, prefix="--bool-true", default=True),
            ToolInput("flag_false", Boolean, position=5, prefix="--bool-false", default=False),
            
            ToolInput("opt_basic", String, position=6, prefix="--opt-basic"),
            ToolInput("opt_default", String, position=7, default="DEFAULT", prefix="--opt-default"),
            ToolInput("opt_optional", String(optional=True), position=8, prefix="--opt-optional"),
        ]

    def outputs(self):
        return [ToolOutput("out", Stdout)]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"


class SecondariesTestTool(CommandTool):
    def tool(self) -> str:
        return "SecondariesTestTool"

    def base_command(self) -> Optional[str | list[str]]:
        return []

    def inputs(self) -> list[ToolInput]:
        return [ToolInput("inp", BamBai, position=4)]

    def arguments(self) -> list[ToolArgument]:
        return [
            ToolArgument("echo hello > out.bam", position=1, shell_quote=False),
            ToolArgument("&& echo there > out.bam.bai", position=2, shell_quote=False),
            ToolArgument("&& echo", position=3, shell_quote=False),
        ]

    def outputs(self):
        return [
            ToolOutput(
                "out", 
                BamBai,
                selector=WildcardSelector("*.bam"),
                secondaries_present_as={".bai": ".bai"},
            )
        ]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"



