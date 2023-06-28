
from typing import Optional

from janis_core import (
    Workflow,
    CommandTool,
    ToolInput,
    ToolOutput,
    InputSelector
)
from janis_core.types import (
    Stdout,
    String,
    File,
)


class Subworkflow2TestWF(Workflow):
    def constructor(self):
        self.input('inFile', File)
        self.input('newFilename', String(optional=True))

        self.step(
            "stp1", 
            GetFilename(inp=self.inFile)
        )
        self.step(
            "stp2", 
            RenameFile(
                inp=self.inFile,
                newFilename=self.newFilename,
            )
        )
        self.step(
            "stp3", 
            SubWF(
                inFile=self.inFile,
                newFilename=self.newFilename,
            )
        )

    def friendly_name(self):
        return "TEST: Subworkflow2TestWF"

    def id(self) -> str:
        return self.__class__.__name__
    

class Subworkflow3TestWF(Workflow):
    def constructor(self):
        self.input('inFile', File)
        self.input('newFilename', String(optional=True))

        self.step(
            "stp2", 
            RenameFile(
                inp=self.inFile,
                newFilename=self.newFilename,
            )
        )
        self.step(
            "stp4", 
            SubWF(
                inFile=self.inFile,
                newFilename=self.stp2.out,
            )
        )

    def friendly_name(self):
        return "TEST: Subworkflow3TestWF"

    def id(self) -> str:
        return self.__class__.__name__
    


class GetFilename(CommandTool):
    def tool(self) -> str:
        return "GetFilename"

    def base_command(self) -> Optional[str | list[str]]:
        return "echo"

    def inputs(self) -> list[ToolInput]:
        return [ToolInput("inp", File, position=1)]

    def outputs(self):
        return [ToolOutput("out", Stdout)]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"
    

class RenameFile(CommandTool):
    def tool(self) -> str:
        return "RenameFile"

    def base_command(self) -> Optional[str | list[str]]:
        return "cat"

    def inputs(self) -> list[ToolInput]:
        return [
            ToolInput("inp", File, position=1),
            ToolInput("newFilename", String(optional=True), prefix='>', position=2)
        ]

    def outputs(self):
        return [ToolOutput("out", File(optional=True), selector=InputSelector('newFilename'))]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"



class SubWF(Workflow):
    def constructor(self):
        self.input('inFile', File)
        self.input('newFilename', String(optional=True))

        self.step(
            "stp1", 
            GetFilename(inp=self.inFile)
        )
        self.step(
            "stp2", 
            RenameFile(
                inp=self.inFile,
                newFilename=self.newFilename,
            )
        )

    def friendly_name(self):
        return "TEST: SubWF"

    def id(self) -> str:
        return self.__class__.__name__