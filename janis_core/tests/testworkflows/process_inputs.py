
from janis_core import Workflow, ToolInput, ToolOutput, CommandTool
from janis_core.types import Array
from janis_core.redefinitions.types import Fastq
from typing import Optional



class ProcessInputsTestWF(Workflow):

    def constructor(self):
        self.input('inFastq', Fastq)
        self.input('inFastqArray', Array(Fastq))

        self.step(
            "stp1",
            StageAsTestTool1(
                fastq_inp1=self.inFastq,
                fastq_inp2=self.inFastq,
            )
        )

        self.step(
            "stp2",
            StageAsTestTool1(
                fastq_inp1=self.stp1.fastq_out1,
                fastq_inp2=self.stp1.fastq_out2,
            )
        )

    def friendly_name(self):
        return "TEST: PlumbingTypeMismatchTestWF"

    def id(self) -> str:
        return self.__class__.__name__



class StageAsTestTool1(CommandTool):
    def tool(self) -> str:
        return "StageAsTestTool"
    
    def base_command(self) -> Optional[str | list[str]]:
        return "echo"
    
    def inputs(self) -> list[ToolInput]:
        return [
            ToolInput("fastq_inp1", Fastq, position=0),
            ToolInput("fastq_inp2", Fastq, position=1),
        ]

    def outputs(self):
        return [
            ToolOutput("fastq_out1", Fastq(), selector='dir1/fastq.fq'),
            ToolOutput("fastq_out2", Fastq(), selector='dir2/fastq.fq'),
        ]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"

