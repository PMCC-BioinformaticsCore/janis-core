

from typing import Optional

from janis_core import Workflow, ToolInput, ToolOutput, CommandTool
from janis_core.types import File, Array, Stdout
from janis_core.redefinitions.types import BamBai, Bam, FastaWithIndexes, FastaDict


class PlumbingTypeMismatchTestWF(Workflow):

    def constructor(self):
        self.input('inFile', File)
        self.input('inFileArray', Array(File))
        self.input('inBamBai', BamBai)
        self.input('inBamBaiArray', Array(BamBai))
        self.input('inFastaWithIndexes', FastaWithIndexes)

        # single/array mismatches
        self.step(
            "single_to_array", 
            ArrayInputTestTool(inp=self.inFile)
        )
        self.step(
            "array_to_single",
            SingleInputTestTool(inp=self.inFileArray)
        )

        # secondary/secondary mismatches
        self.step(
            "bambai_to_bam",
            BamInputTestTool(inp=self.inBamBai)
        )
        self.step(
            "bambai_to_bam_array",
            BamArrayInputTestTool(inp=self.inBamBai)
        )

        # secondary/secondary array mismatches
        self.step(
            "fastawithindexes_to_fastadict",
            FastaDictInputTestTool(inp=self.inFastaWithIndexes)
        )
        self.step(
            "fastawithindexes_to_fastadict_array",
            FastaDictArrayInputTestTool(inp=self.inFastaWithIndexes)
        )

        # secondary array mismatches
        self.step(
            "secondary_array_to_secondary", 
            SecondaryInputTestTool(inp=self.inBamBaiArray)
        )
        self.step(
            "secondary_array_to_secondary_array", 
            SecondaryArrayInputTestTool(inp=self.inBamBaiArray)
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


class FastaDictInputTestTool(BaseTestTool):
    
    def tool(self) -> str:
        return "FastaDictInputTestTool"
    
    def inputs(self) -> list[ToolInput]:
        return [ToolInput("inp", FastaDict(), position=0)]


class FastaDictArrayInputTestTool(BaseTestTool):
    
    def tool(self) -> str:
        return "FastaDictArrayInputTestTool"
    
    def inputs(self) -> list[ToolInput]:
        return [ToolInput("inp", Array(FastaDict()), position=0)]


class BamInputTestTool(BaseTestTool):
    
    def tool(self) -> str:
        return "BamInputTestTool"
    
    def inputs(self) -> list[ToolInput]:
        return [ToolInput("inp", Bam(), position=0)]


class BamArrayInputTestTool(BaseTestTool):
    
    def tool(self) -> str:
        return "BamInputTestTool"
    
    def inputs(self) -> list[ToolInput]:
        return [ToolInput("inp", Array(Bam()), position=0)]


class SingleInputTestTool(BaseTestTool):
    
    def tool(self) -> str:
        return "SingleInputTestTool"
    
    def inputs(self) -> list[ToolInput]:
        return [ToolInput("inp", File(), position=0)]


class ArrayInputTestTool(BaseTestTool):
    
    def tool(self) -> str:
        return "ArrayInputTestTool"
    
    def inputs(self) -> list[ToolInput]:
        return [ToolInput("inp", Array(File()), position=0)]


class SecondaryInputTestTool(BaseTestTool):
    
    def tool(self) -> str:
        return "SecondaryInputTestTool"
    
    def inputs(self) -> list[ToolInput]:
        return [ToolInput("inp", BamBai(), position=0)]


class SecondaryArrayInputTestTool(BaseTestTool):
    
    def tool(self) -> str:
        return "SecondaryArrayInputTestTool"
    
    def inputs(self) -> list[ToolInput]:
        return [ToolInput("inp", Array(BamBai()), position=0)]


