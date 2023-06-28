

from typing import Optional
from janis_bioinformatics.data_types.bam import BamBai
from janis_core.types import Array, Stdout
from janis_core import (
    Workflow,
    CommandTool,
    ToolInput,
    ToolOutput,
)


class SecondaryArrayToolInput(Workflow):
    def id(self) -> str:
        return "SecondaryArrayToolInputTestWF"

    def friendly_name(self):
        return "TEST: SecondaryArrayToolInputTestWF"

    def constructor(self):
        self.input('inBamBaiArr', Array(BamBai()))

        self.step(
            "stp1", 
            TestTool(
                inBamBaiArr=self.inBamBaiArr
            ), 
        )
        self.step(
            "stp2", 
            TestTool(
                inBamBai=self.inBamBaiArr
            ), 
            scatter="inBamBai"
        )

        self.output("out1", source=self.stp1.out)
        self.output("out2", source=self.stp2.out)



# TOOLS

class TestTool(CommandTool):
    def tool(self) -> str:
        return "TestTool"

    def base_command(self) -> Optional[str | list[str]]:
        return ['echo']

    def inputs(self) -> list[ToolInput]:
        return [
            ToolInput(
                "inp", 
                Array(BamBai, optional=True), 
                position=1
            ),
            ToolInput(
                "inp", 
                BamBai(optional=True), 
                position=1
            ),
        ]

    def outputs(self):
        return [
            ToolOutput(
                "out", 
                Stdout(),
            )
        ]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"

