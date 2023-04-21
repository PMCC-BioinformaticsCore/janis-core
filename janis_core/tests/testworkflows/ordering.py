from typing import Optional

from janis_core import (
    Workflow,
    CommandTool,
    ToolInput,
    ToolOutput,
    InputSelector
)

from janis_core.types import (
    String,
    File,
    Array,
    Int,
)

from janis_bioinformatics.data_types.fastq import Fastq



class OrderingTestWF(Workflow):

    def constructor(self):
        self.input('inFile', File())
        self.input('inIntArray', Array(Int()))
        self.input('inStr', String())
        self.input('inInt', Int())
        self.input('inFastq', Fastq())
        self.input('inFastqArray', Array(Fastq()))

        # from workflow inputs
        self.step(
            "stp1", 
            MultiTypeTestTool(
                inFastq=self.inFastq,
                inFastqArray=self.inFastqArray,
                inFile=self.inFile,
                inStr=self.inStr,
                inInt=self.inInt,
                inIntArray=self.inIntArray,
            )
        )
        self.step(
            "stp2", 
            MultiTypeTestWF(
                inFastq=self.inFastq,
                inFastqArray=self.inFastqArray,
                inFile=self.inFile,
                inStr=self.inStr,
                inInt=self.inInt,
                inIntArray=self.inIntArray,
            )
        )

        # from process outputs
        self.step(
            "stp3", 
            MultiTypeTestTool(
                inFastq=self.stp1.outFastq,
                inFastqArray=self.stp1.outFastqArray,
                inFile=self.stp1.outFile,
                inStr=self.stp1.outStr,
                inInt=self.stp1.outInt,
                inIntArray=self.stp1.outIntArray,
            )
        )
        self.step(
            "stp4", 
            MultiTypeTestWF(
                inFastq=self.stp1.outFastq,
                inFastqArray=self.stp1.outFastqArray,
                inFile=self.stp1.outFile,
                inStr=self.stp1.outStr,
                inInt=self.stp1.outInt,
                inIntArray=self.stp1.outIntArray,
            )
        )
        
        # from subworkflow outputs
        self.step(
            "stp5", 
            MultiTypeTestTool(
                inFastq=self.stp2.outFastq,
                inFastqArray=self.stp2.outFastqArray,
                inFile=self.stp2.outFile,
                inStr=self.stp2.outStr,
                inInt=self.stp2.outInt,
                inIntArray=self.stp2.outIntArray,
            )
        )
        self.step(
            "stp6", 
            MultiTypeTestWF(
                inFastq=self.stp2.outFastq,
                inFastqArray=self.stp2.outFastqArray,
                inFile=self.stp2.outFile,
                inStr=self.stp2.outStr,
                inInt=self.stp2.outInt,
                inIntArray=self.stp2.outIntArray,
            )
        )

    def friendly_name(self):
        return "TEST: OrderingTestWF"

    def id(self) -> str:
        return self.__class__.__name__



class MultiTypeTestTool(CommandTool):
    def tool(self) -> str:
        return "MultiTypeTestTool"

    def base_command(self) -> Optional[str | list[str]]:
        return "echo"

    def inputs(self) -> list[ToolInput]:
        return [
            ToolInput("inStr", String(), position=4),
            ToolInput("inFastq", Fastq(), position=1),
            ToolInput("inIntArray", Array(Int()), position=6),
            ToolInput("inFastqArray", Array(Fastq()), position=2),
            ToolInput("inInt", Int(), position=5),
            ToolInput("inFile", File(), position=3),
        ]

    def outputs(self):
        return [
            ToolOutput("outFastq", Fastq(), selector=InputSelector('inFastq')),
            ToolOutput("outFastqArray", Array(Fastq()), selector=InputSelector('inFastqArray')),
            ToolOutput("outFile", File(), selector=InputSelector('inFile')),
            ToolOutput("outStr", String(), selector=InputSelector('inStr')),
            ToolOutput("outInt", Int(), selector=InputSelector('inInt')),
            ToolOutput("outIntArray", Array(Int()), selector=InputSelector('inIntArray')),
        ]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"



# all WildcardSelector use cases
class MultiTypeTestWF(Workflow):

    def constructor(self):
        self.input('inStr', String)
        self.input('inIntArray', Array(Int))
        self.input('inFastq', Fastq())
        self.input('inInt', Int)
        self.input('inFastqArray', Array(Fastq()))
        self.input('inFile', File)

        self.output('outFastq', Fastq(), source=self.inFastq)
        self.output('outFastqArray', Array(Fastq()), source=self.inFastqArray)
        self.output('outFile', File, source=self.inFile)
        self.output('outStr', String, source=self.inStr)
        self.output('outInt', Int, source=self.inInt)
        self.output('outIntArray', Array(Int), source=self.inIntArray)

    def friendly_name(self):
        return "TEST: MultiTypeTestWF"

    def id(self) -> str:
        return self.__class__.__name__




