

from typing import Optional

from janis_core.types import File, Array, Stdout, String
from janis_core import (
    Workflow, 
    ToolInput, 
    ToolOutput, 
    CommandTool,
)

from janis_core import (
    If, 
    IsDefined,
    FirstOperator
)


class EntityTraceTestWF(Workflow):

    def constructor(self):
        # self.input("inString", String(optional=True), value="")
        self.input("inFile", File())
        self.input("inFileArray", Array(File()))
        self.input("inString", String())
        self.input("inStringOpt", String(optional=True))
        self.input("inStringOptBackup", String(optional=True))
        self.input("inStringArray", Array(String()))

        # File -> File workflow input
        self.step(
            "stp1",
            CatTestTool(inp=self.inFile)
        )
        
        # File -> File step connection
        self.step(
            "stp2",
            CatTestTool(inp=self.stp1.out)
        )
        
        # If & IsDefined operators
        self.step(
            "stp3",
            EchoTestTool(inp=If(IsDefined(self.inStringOpt), self.inStringOpt, self.inStringOptBackup)),
        )
        
        # FirstOperator
        derived_string = FirstOperator(
            [
                self.inStringOpt,
                self.step(
                    "get_string",
                    EchoTestTool(inp="Some default value"),
                    when=self.inStringOpt.is_null(),
                ).out,
            ]
        )
        self.step(
            "stp4",
            EchoTestTool(inp=derived_string),
        )

        # IndexOperator
        self.step(
            "stp5", 
            EchoTestTool(inp=self.inStringArray[0])
        )

        self.step(
            "stp6",
            CatTestTool(inp=self.inFileArray),
            scatter='inp'
        )
        
        self.step(
            "stp7",
            CatTestTool(inp=self.stp6.out),
        )


    def friendly_name(self):
        return "TEST: EntityTraceTestWF"

    def id(self) -> str:
        return self.__class__.__name__


class EchoTestTool(CommandTool):
    def base_command(self) -> Optional[str | list[str]]:
        return "echo"
    
    def tool(self) -> str:
        return "EchoTestTool"

    def inputs(self) -> list[ToolInput]:
        return [ToolInput("inp", String, position=0)]

    def outputs(self):
        return [ToolOutput("out", Stdout())]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"


class CatTestTool(CommandTool):
    def base_command(self) -> Optional[str | list[str]]:
        return "echo"
    
    def tool(self) -> str:
        return "CatTestTool"

    def inputs(self) -> list[ToolInput]:
        return [ToolInput("inp", File, position=0)]

    def outputs(self):
        return [ToolOutput("out", Stdout())]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"


