

from typing import Optional

from janis_core import (
    Workflow,
    CommandTool,
    ToolInput,
    ToolOutput,
    WildcardSelector,
    InputSelector
)

from janis_core.types import (
    String,
    File,
    Array,
    Stdout
)

from janis_bioinformatics.data_types.bam import BamBai
from janis_bioinformatics.data_types.vcf import VcfTabix
from janis_bioinformatics.data_types.fasta import FastaWithIndexes


class NamingTestWF(Workflow):

    def constructor(self):
        self.input('processInput', File())
        self.input('paramInput', String())
        self.input('secondary', BamBai())
        self.input('processInputArray', Array(File()))
        self.input('paramInputArray', Array(String()))
        self.input('normal_sample', Array(BamBai()))
        self.input('tumor_sample', Array(BamBai()))
        self.input('reference', FastaWithIndexes())
        self.input('panelOfNormals', VcfTabix())
        self.input('germlineResource', VcfTabix())

        self.step(
            "stp1", 
            NamingTestTool( 
                processInput=self.processInput,
                paramInput=self.paramInput,
                secondary=self.secondary,
                processInputArray=self.processInputArray,
                paramInputArray=self.paramInputArray,
                secondaryArray=self.normal_sample,
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
                secondaryArray=self.normal_sample,
            )
        )
        
        self.step(
            "stp3", 
            HardNamingTestTool(
                reference=self.reference,
                panelOfNormals=self.panelOfNormals,
                germlineResource=self.germlineResource,
                normal_samples=self.normal_sample,
                tumor_samples=self.tumor_sample,
            )
        )

    def friendly_name(self):
        return "TEST: NamingTestWF"

    def id(self) -> str:
        return self.__class__.__name__



class HardNamingTestTool(CommandTool):
    def tool(self) -> str:
        return "HardNamingTestTool"

    def base_command(self) -> Optional[str | list[str]]:
        return "echo"

    def inputs(self) -> list[ToolInput]:
        return [
            ToolInput('reference', FastaWithIndexes(), prefix='--reference', position=1),
            ToolInput('panelOfNormals', VcfTabix(), prefix='--normals', position=2),
            ToolInput('germlineResource', VcfTabix(), prefix='--germline', position=3),
            ToolInput('normal_samples', Array(BamBai()), prefix='--normal_samples', position=4),
            ToolInput('tumor_samples', Array(BamBai()), prefix='--tumor_samples', position=5),
        ]

    def outputs(self):
        return [ToolOutput("out", Stdout())]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"




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
            ToolOutput("outProcessInput", File(), selector=InputSelector('processInput') + '.fastq'),
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




