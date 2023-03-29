

from janis_core import Workflow

from janis_core.types import (
    File,
    Array,
    Stdout, 
    Filename
)
from typing import Optional
from janis_core import CommandTool, ToolInput, ToolOutput

from janis_bioinformatics.data_types.fastq import FastqGzPair
from janis_bioinformatics.data_types.bam import BamBai


# TODO add filename test

class OptionalTestWF(Workflow):

    def constructor(self):
        self.input('inSecondaryArray', Array(BamBai()))
        self.input('inSecondaryArrayOpt', Array(BamBai(), optional=True))
        
        self.input('inSecondary', BamBai())
        self.input('inSecondaryOpt', BamBai(optional=True))

        self.input('inFilePair', FastqGzPair())
        self.input('inFilePairOpt', FastqGzPair(optional=True))

        self.input('inFileArray', Array(File()))
        self.input('inFileArrayOpt', Array(File(), optional=True))

        self.input('inFile', File())
        self.input('inFileOpt', File(optional=True))

        self.step(
            "stp1", 
            OptionalTestTool(
                inSecondaryArrayOpt=self.inSecondaryArrayOpt,
                inSecondaryOpt=self.inSecondaryOpt,
                inFilePairOpt=self.inFilePairOpt,
                inFileArrayOpt=self.inFileArrayOpt,
                inFileOpt=self.inFileOpt,
            )
        )
        self.step(
            "stp2", 
            MandatoryTestTool(
                inSecondaryArray=self.inSecondaryArray,
                inSecondary=self.inSecondary,
                inFilePair=self.inFilePair,
                inFileArray=self.inFileArray,
                inFile=self.inFile,
            )
        )

    def friendly_name(self):
        return "TEST: OptionalTestWF"

    def id(self) -> str:
        return self.__class__.__name__



class OptionalTestTool(CommandTool):
    
    def tool(self) -> str:
        return "OptionalTestTool"
        
    def base_command(self) -> Optional[str | list[str]]:
        return "echo"

    def inputs(self) -> list[ToolInput]:
        return [
            ToolInput('inSecondaryArrayOpt', Array(BamBai(), optional=True), position=1),
            ToolInput('inSecondaryOpt', BamBai(optional=True), position=2),
            ToolInput('inFilePairOpt', FastqGzPair(optional=True), position=3),
            ToolInput('inFileArrayOpt', Array(File(), optional=True), position=4),
            ToolInput('inFileOpt', File(optional=True), position=5),
        ]

    def outputs(self):
        return [ToolOutput('out', Stdout())]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"
    

class MandatoryTestTool(CommandTool):
    
    def tool(self) -> str:
        return "MandatoryTestTool"
        
    def base_command(self) -> Optional[str | list[str]]:
        return "echo"

    def inputs(self) -> list[ToolInput]:
        return [
            ToolInput('inSecondaryArray', Array(BamBai()), position=1),
            ToolInput('inSecondary', BamBai(), position=2),
            ToolInput('inFilePair', FastqGzPair(), position=3),
            ToolInput('inFileArray', Array(File()), position=4),
            ToolInput('inFile', File(), position=5),
        ]

    def outputs(self):
        return [ToolOutput('out', Stdout())]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"

