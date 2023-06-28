

from typing import Optional

from janis_core import Workflow
from janis_core import CommandTool, ToolInput, ToolOutput
from janis_core.types import (
    File,
    Array,
    Stdout, 
    Int,
    String
)

from janis_bioinformatics.data_types.fastq import FastqGzPair
from janis_bioinformatics.data_types.bam import BamBai


class MandatoryInputTypesTestWF(Workflow):

    def constructor(self):
        self.input('inSecondaryArray', Array(BamBai()))
        self.input('inSecondary', BamBai())
        self.input('inFilePairArray', Array(FastqGzPair()))
        self.input('inFilePair', FastqGzPair())
        self.input('inFileArray', Array(File()))
        self.input('inFile', File())
        self.input('inIntArray', Array(Int()))
        self.input('inInt', Int())
        self.input('inStrArray', Array(String()))
        self.input('inStr', String())

        self.step(
            "stp1", 
            MandatoryInputTypesTestTool(
                inSecondaryArray=self.inSecondaryArray,
                inSecondary=self.inSecondary,
                inFilePairArray=self.inFilePairArray,
                inFilePair=self.inFilePair,
                inFileArray=self.inFileArray,
                inFile=self.inFile,
                inIntArray=self.inIntArray,
                inInt=self.inInt,
                inStrArray=self.inStrArray,
                inStr=self.inStr
            )
        )
        self.step(
            "stp2", 
            MandatoryInputTypesTestTool(
                inSecondaryArray=self.inSecondaryArray,
                inSecondary=self.inSecondary,
                inFilePairArray=self.inFilePairArray,
                inFilePair=self.inFilePair,
                inFileArray=self.inFileArray,
                inFile=self.inFile,
                inIntArray=[1, 2, 3],
                inInt=10,
                inStrArray=['hi', 'there'],
                inStr='friend'
            )
        )

    def friendly_name(self):
        return "TEST: MandatoryInputTypesTestWF"

    def id(self) -> str:
        return self.__class__.__name__


class MandatoryInputTypesTestTool(CommandTool):
    
    def tool(self) -> str:
        return "MandatoryInputTypesTestTool"
        
    def base_command(self) -> Optional[str | list[str]]:
        return "echo"

    def inputs(self) -> list[ToolInput]:
        return [
            ToolInput('inSecondaryArray', Array(BamBai()), position=1),
            ToolInput('inSecondary', BamBai(), position=2),
            ToolInput('inFilePairArray', Array(FastqGzPair()), position=3),
            ToolInput('inFilePair', FastqGzPair(), position=4),
            ToolInput('inFileArray', Array(File()), position=5),
            ToolInput('inFile', File(), position=6),
            ToolInput('inIntArray', Array(Int()), position=7),
            ToolInput('inInt', Int(), position=8),
            ToolInput('inStrArray', Array(String()), position=9),
            ToolInput('inStr', String(), position=10),
        ]

    def outputs(self):
        return [ToolOutput('out', Stdout())]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"

