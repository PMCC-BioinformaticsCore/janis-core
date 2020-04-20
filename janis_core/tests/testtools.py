from typing import Dict, List, Any, Optional, Union
from janis_core import (
    ToolOutput,
    ToolInput,
    String,
    CommandTool,
    Stdout,
    InputSelector,
    Array,
    File,
    WildcardSelector,
    StringFormatter,
    ToolArgument,
    InputDocumentation,
    InputQualityType,
)


class EchoTestTool(CommandTool):
    def tool(self) -> str:
        return "EchoTestTool"

    def base_command(self) -> Optional[Union[str, List[str]]]:
        return "echo"

    def inputs(self) -> List[ToolInput]:
        return [ToolInput("inp", str, position=0)]

    def outputs(self):
        return [ToolOutput("out", Stdout)]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"


class CatTestTool(CommandTool):
    def tool(self) -> str:
        return "CatTestTool"

    def base_command(self) -> Optional[Union[str, List[str]]]:
        return "cat"

    def inputs(self) -> List[ToolInput]:
        return [ToolInput("inp", File, position=0)]

    def outputs(self):
        return [ToolOutput("out", Stdout)]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"


class SingleTestTool(CommandTool):
    @staticmethod
    def tool():
        return "TestStepTool"

    @staticmethod
    def base_command():
        return "echo"

    def inputs(self):
        return [
            ToolInput("inputs", String(), position=0),
            ToolInput("input2", String(optional=True), position=1),
            ToolInput("input3", String(optional=True), position=2),
            ToolInput("input4", String(optional=True), position=3),
        ]

    def friendly_name(self):
        return None

    def outputs(self):
        return [ToolOutput("out", String(), glob=WildcardSelector("*"))]

    @staticmethod
    def container():
        return None

    @staticmethod
    def version():
        return None


class ArrayTestTool(CommandTool):
    @staticmethod
    def tool():
        return "ArrayStepTool"

    def friendly_name(self):
        return None

    @staticmethod
    def base_command():
        return "echo"

    def inputs(self):
        return [ToolInput("inputs", Array(String()))]

    def outputs(self):
        return [ToolOutput("outs", Array(File()), glob=WildcardSelector("*"))]

    @staticmethod
    def container():
        return None

    @staticmethod
    def version():
        return None


class TestTool(CommandTool):
    @staticmethod
    def tool():
        return "TestTranslationtool"

    @staticmethod
    def base_command():
        return "echo"

    def inputs(self) -> List[ToolInput]:
        return [ToolInput("testtool", String())]

    def arguments(self) -> List[ToolArgument]:
        return [ToolArgument(StringFormatter('test:\\t:escaped:\\n:characters"'))]

    def outputs(self) -> List[ToolOutput]:
        return [ToolOutput("std", Stdout())]

    def cpus(self, hints: Dict[str, Any]):
        return 2

    def memory(self, hints: Dict[str, Any]):
        return 2

    def friendly_name(self) -> str:
        return "Tool for testing translation"

    @staticmethod
    def container():
        return "ubuntu:latest"

    @staticmethod
    def version():
        return None

    def env_vars(self):
        return {"test1": InputSelector("testtool")}


class TestToolWithSecondaryOutput(TestTool):
    def outputs(self):
        return [
            ToolOutput(
                "out", TestTypeWithSecondary(), glob=InputSelector("testtool") + "/out"
            )
        ]


class TestToolWithSecondaryInput(CatTestTool):
    def inputs(self) -> List[ToolInput]:
        return [ToolInput("inp", TestTypeWithSecondary, position=0)]


class TestTypeWithSecondary(File):
    @staticmethod
    def secondary_files():
        return ["^.txt"]


class TestTypeWithNonEscapedSecondary(File):
    @staticmethod
    def secondary_files():
        return [".txt"]


class TestInputQualityTool(CommandTool):
    def tool(self) -> str:
        return "TESTONLY_inputQualityTool"

    def base_command(self) -> Optional[Union[str, List[str]]]:
        return "echo"

    def inputs(self) -> List[ToolInput]:
        return [
            ToolInput(
                "user", str, doc=InputDocumentation(None, quality=InputQualityType.user)
            ),
            ToolInput(
                "static",
                str,
                doc=InputDocumentation(None, quality=InputQualityType.static),
            ),
            ToolInput(
                "configuration",
                str,
                doc=InputDocumentation(None, quality=InputQualityType.configuration),
            ),
            ToolInput("none", str, doc=InputDocumentation(None, quality=None)),
        ]

    def outputs(self) -> List[ToolOutput]:
        return [ToolOutput("out", Stdout)]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"
