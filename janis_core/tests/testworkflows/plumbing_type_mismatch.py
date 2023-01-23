

from typing import Optional

from janis_core import Workflow, ToolInput, ToolOutput, CommandTool
from janis_core.types import File, Array, Stdout
from janis_bioinformatics.data_types.bam import BamBai


class PlumbingTypeMismatchTestWF(Workflow):

    def constructor(self):
        self.input('inFile', File)
        self.input('inFileArray', Array(File))
        self.input('inSecondary', BamBai)
        self.input('inSecondaryArray', Array(BamBai))

        self.step(
            "array_to_single",
            SingleInputTestTool(inp=self.inFileArray)
        )
        self.step(
            "single_to_array", 
            ArrayInputTestTool(inp=self.inFile)
        )
        self.step(
            "secondary_array_to_secondary", 
            SecondaryInputTestTool(inp=self.inSecondaryArray)
        )
        self.step(
            "secondary_array_to_secondary_array", 
            SecondaryArrayInputTestTool(inp=self.inSecondaryArray)
        )

    def friendly_name(self):
        return "TEST: PlumbingTypeMismatchTestWF"

    def id(self) -> str:
        return self.__class__.__name__



class BaseTestTool(CommandTool):
    def base_command(self) -> Optional[str | list[str]]:
        return "echo"

    def outputs(self):
        return [ToolOutput("out", Stdout())]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"


class SingleInputTestTool(BaseTestTool):
    
    def tool(self) -> str:
        return "SingleInputTestTool"
    
    def inputs(self) -> list[ToolInput]:
        return [ToolInput("inp", File, position=0)]


class ArrayInputTestTool(BaseTestTool):
    
    def tool(self) -> str:
        return "ArrayInputTestTool"
    
    def inputs(self) -> list[ToolInput]:
        return [ToolInput("inp", Array(File), position=0)]


class SecondaryInputTestTool(BaseTestTool):
    
    def tool(self) -> str:
        return "SecondaryInputTestTool"
    
    def inputs(self) -> list[ToolInput]:
        return [ToolInput("inp", BamBai, position=0)]


class SecondaryArrayInputTestTool(BaseTestTool):
    
    def tool(self) -> str:
        return "SecondaryArrayInputTestTool"
    
    def inputs(self) -> list[ToolInput]:
        return [ToolInput("inp", Array(BamBai), position=0)]


