

from typing import Optional 

from janis_bioinformatics.data_types.bam import BamBai
from janis_core import (
    ToolOutput,
    ToolInput,
    CommandTool,
    WildcardSelector,
    ToolArgument
)

from janis_core.types import (
    File,
    String,
    Int,
    Stdout,
    Array,
    Boolean
)


class ArrayFileTestTool(CommandTool):
    def tool(self) -> str:
        return "InFileArrayTestTool"

    def base_command(self) -> Optional[str | list[str]]:
        return "echo"

    def inputs(self) -> list[ToolInput]:
        return [ToolInput("ins", Array(File), position=0)]

    def arguments(self) -> Optional[list[ToolArgument]]:
        return [
            ToolArgument("&& echo hi > hi.txt", position=1, shell_quote=False),
            ToolArgument("&& echo there > there.txt", position=2, shell_quote=False),
            ToolArgument("&& echo friend > friend.txt", position=3, shell_quote=False),
        ]

    def outputs(self):
        return [ToolOutput("outs", Array(File()), selector=WildcardSelector("*.txt"))]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"


class ArrayStringTestTool(CommandTool):
    def tool(self):
        return "InStringArrayTestTool"

    def friendly_name(self):
        return None

    def base_command(self):
        return "echo"

    def inputs(self):
        return [ToolInput("ins", Array(String()), position=1)]

    def outputs(self):
        return [ToolOutput("outs", Array(File()), glob=WildcardSelector("*"))]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"


class ArrayIntTestTool(CommandTool):
    def tool(self):
        return "InIntArrayTestTool"

    def friendly_name(self):
        return None

    def base_command(self):
        return "echo"

    def inputs(self):
        return [ToolInput("ins", Array(Int()), position=1)]

    def outputs(self):
        return [ToolOutput("outs", Array(File()), glob=WildcardSelector("*"))]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"


class ArrayComponentsTestTool(CommandTool):
    def tool(self) -> str:
        return "ArrayComponentTypeTestTool"

    def base_command(self) -> Optional[str | list[str]]:
        return "echo"

    def inputs(self) -> list[ToolInput]:
        return [
            ToolInput("pos_basic", Array(File), position=1),
            ToolInput("pos_default", Array(String), default=["hi", "there", "friend"], position=2),
            ToolInput("pos_optional", Array(String, optional=True), position=3),

            ToolInput("flag_true", Array(Boolean), position=4, prefix="--bool-true", default=[True]),
            ToolInput("flag_false", Array(Boolean), position=5, prefix="--bool-false", default=[True]),
            
            ToolInput("opt_basic", Array(String), position=6, prefix="--opt-basic=", separate_value_from_prefix=False),
            ToolInput("opt_default", Array(String), position=7, default=["hi", "there", "friend"], prefix="--opt-default", prefix_applies_to_all_elements=True),
            ToolInput("opt_optional", Array(String, optional=True), position=8, prefix="--opt-optional", separator=","),
        ]

    def outputs(self):
        return [ToolOutput("out", Stdout)]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"


# TODO FAILS secondaries_present_as
class ArraySecondariesTestTool(CommandTool):
    def tool(self) -> str:
        return "ArraySecondariesTestTool"

    def base_command(self) -> Optional[str | list[str]]:
        return []

    def inputs(self) -> list[ToolInput]:
        return [ToolInput("inp", Array(BamBai), position=6)]

    def arguments(self) -> list[ToolArgument]:
        return [
            ToolArgument("echo 1 > out1.bam", position=1, shell_quote=False),
            ToolArgument("&& echo 1 > out1.bam.bai", position=2, shell_quote=False),
            ToolArgument("&& echo 2 > out2.bam", position=3, shell_quote=False),
            ToolArgument("&& echo 2 > out2.bam.bai", position=4, shell_quote=False),
            ToolArgument("&& echo", position=5, shell_quote=False),
        ]

    def outputs(self):
        return [
            ToolOutput("outStdout", Stdout),
            ToolOutput(
                "outArray", 
                Array(BamBai),
                selector=WildcardSelector("*.bam"),
                secondaries_present_as={".bai": ".bai"},
            )
        ]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"








### this stuff is old. leave. ###


class ArrayStepTool(CommandTool):
    def tool(self):
        return "ArrayStepTool"

    def friendly_name(self):
        return None

    def base_command(self):
        return "echo"

    def inputs(self):
        return [ToolInput("inps", Array(String()), position=1)]

    def outputs(self):
        return [ToolOutput("outs", Array(File()), glob=WildcardSelector("*"))]

    def container(self):
        return None

    def version(self):
        return None







