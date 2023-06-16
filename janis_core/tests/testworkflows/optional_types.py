

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
from janis_core.redefinitions.types import BamBai, FastqGzPair


class OptionalInputTypesTestWF(Workflow):

    def constructor(self):
        self.input('inSecondaryArray', Array(BamBai(), optional=True))
        self.input('inSecondary', BamBai(optional=True))
        self.input('inFilePairArray', Array(FastqGzPair(), optional=True))
        self.input('inFilePair', FastqGzPair(optional=True))
        self.input('inFileArray', Array(File(), optional=True))
        self.input('inFile', File(optional=True))
        self.input('inIntArray', Array(Int(), optional=True))
        self.input('inInt', Int(optional=True))
        self.input('inStrArray', Array(String(), optional=True))
        self.input('inStr', String(optional=True))

        self.step(
            "stp1", 
            OptionalInputTypesTestTool(
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
            OptionalInputTypesTestTool(
                inSecondaryArray=None,
                inSecondary=None,
                inFilePairArray=None,
                inFilePair=None,
                inFileArray=None,
                inFile=None,
                inIntArray=None,
                inInt=None,
                inStrArray=None,
                inStr=None,
            )
        )

    def friendly_name(self):
        return "TEST: OptionalInputTypesTestWF"

    def id(self) -> str:
        return self.__class__.__name__


class OptionalInputTypesTestTool(CommandTool):
    
    def tool(self) -> str:
        return "OptionalInputTypesTestTool"
        
    def base_command(self) -> Optional[str | list[str]]:
        return "echo"

    def inputs(self) -> list[ToolInput]:
        return [
            ToolInput('inSecondaryArray', Array(BamBai(), optional=True), position=1),
            ToolInput('inSecondary', BamBai(optional=True), position=2),
            ToolInput('inFilePairArray', Array(FastqGzPair(), optional=True), position=3),
            ToolInput('inFilePair', FastqGzPair(optional=True), position=4),
            ToolInput('inFileArray', Array(File(), optional=True), position=5),
            ToolInput('inFile', File(optional=True), position=6),
            ToolInput('inIntArray', Array(Int(), optional=True), position=7),
            ToolInput('inInt', Int(optional=True), position=8),
            ToolInput('inStrArray', Array(String(), optional=True), position=9),
            ToolInput('inStr', String(optional=True), position=10),
        ]

    def outputs(self):
        return [ToolOutput('out', Stdout())]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"

