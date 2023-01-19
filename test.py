
from typing import Optional

from janis_core import (
    Workflow,
    CommandTool,
    ToolInput,
    ToolOutput,
)

from janis_core.types import (
    String,
    Array,
    Stdout
)


class TestWF(Workflow):

    def constructor(self):
        self.input('myStringArray', Array(String()))

        self.step(
            "stp1", 
            TestTool( 
                myStringArray=self.myStringArray,
            )
        )

    def friendly_name(self):
        return "TEST: TestWF"

    def id(self) -> str:
        return self.__class__.__name__


class TestTool(CommandTool):
    def tool(self) -> str:
        return "TestTool"

    def base_command(self) -> Optional[str | list[str]]:
        return "echo"

    def inputs(self) -> list[ToolInput]:
        return [
            ToolInput('myStringArray', Array(String()), prefix='--myStringArray', position=1),
        ]

    def outputs(self):
        return [
            ToolOutput("stdout", Stdout())
        ]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"

