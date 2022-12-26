

from typing import Optional

from janis_core import (
    Workflow,
    CommandTool,
    ToolInput,
    ToolOutput,
    WildcardSelector
)

from janis_core.types import (
    String,
    File,
    Array
)

from janis_bioinformatics.data_types.bam import BamBai


class NamingTestWF(Workflow):

    def constructor(self):
        self.input('processInput', File())
        self.input('paramInput', String())
        self.input('secondary', BamBai())
        self.input('processInputArray', Array(File()))
        self.input('paramInputArray', Array(String()))
        self.input('secondaryArray', Array(BamBai()))

        self.step(
            "stp1", 
            NamingTestTool( 
                processInput=self.processInput,
                paramInput=self.paramInput,
                secondary=self.secondary,
                processInputArray=self.processInputArray,
                paramInputArray=self.paramInputArray,
                secondaryArray=self.secondaryArray,
            )
        )
        self.step(
            "stp2", 
            NamingTestTool( 
                processInput=self.stp1.outProcessInput,
                paramInput=self.stp1.outParamInput,
                secondary=self.stp1.outSecondary,
                processInputArray=self.stp1.outProcessInputArray,
                paramInputArray=self.stp1.outParamInputArray,
                secondaryArray=self.secondaryArray,
            )
        )

    def friendly_name(self):
        return "TEST: NamingTestWF"

    def id(self) -> str:
        return self.__class__.__name__


class NamingTestTool(CommandTool):
    def tool(self) -> str:
        return "NamingTestTool"

    def base_command(self) -> Optional[str | list[str]]:
        return "echo"

    def inputs(self) -> list[ToolInput]:
        return [
            ToolInput('processInput', File(), prefix='--processInput', position=1),
            ToolInput('paramInput', String(), prefix='--paramInput', position=2),
            ToolInput('secondary', BamBai(), prefix='--secondary', position=3),
            ToolInput('processInputArray', Array(File()), prefix='--processInputArray', position=4),
            ToolInput('paramInputArray', Array(String()), prefix='--paramInputArray', position=5),
            ToolInput('secondaryArray', Array(BamBai()), prefix='--secondaryArray', position=6),
        ]

    def outputs(self):
        return [
            ToolOutput("outProcessInput", File(), selector=WildcardSelector('process_input.fastq')),
            ToolOutput("outParamInput", File(), selector=WildcardSelector('param_input.txt')),
            ToolOutput(
                "outSecondary", 
                BamBai(), 
                selector=WildcardSelector("*.bam"),
                secondaries_present_as={".bai": ".bai"},
            ),
            ToolOutput("outProcessInputArray", Array(File()), selector=WildcardSelector('process_input_arr*')),
            ToolOutput("outParamInputArray", Array(File()), selector=WildcardSelector('param_input_arr*')),
            # ToolOutput("outSecondaryArray", Array(BamBai()), selector=[WildcardSelector('*.bam'), WildcardSelector('*.bai')]),
        ]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"




