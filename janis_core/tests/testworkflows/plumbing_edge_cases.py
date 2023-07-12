

from typing import Any, Optional

from janis_core import Workflow, PythonTool, CommandTool, TOutput, ToolInput, ToolOutput, File, Stdout
from janis_core.types import Array, String


class PlumbingEdgeCaseTestWF(Workflow):

    def constructor(self):
        self.input('inString', String())

        self.step(
            "stp1",
            PythontoolArrayStringTestTool(inp=self.inString)
        )
        self.step(
            "stp2",
            ArrayStringTestTool(inp=self.stp1.out)
        )

    def friendly_name(self):
        return "TEST: PlumbingEdgeCaseTestWF"

    def id(self) -> str:
        return self.__class__.__name__



class PythontoolArrayStringTestTool(PythonTool):
    @staticmethod
    def code_block(inp: String) -> dict[str, Any]:
        myarray = []
        return {
            'out': myarray
        }

    def tool(self) -> str:
        return "PythontoolArrayStringTestTool"
    
    def outputs(self):
        return [TOutput("out", Array(str))]


class ArrayStringTestTool(CommandTool):
    def base_command(self) -> Optional[str | list[str]]:
        return "echo"
   
    def tool(self) -> str:
        return "ArrayStringTestTool"
    
    def inputs(self) -> list[ToolInput]:
        return [ToolInput("inp", Array(String()), position=0)]

    def outputs(self):
        return [ToolOutput("out", Stdout())]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"


    
