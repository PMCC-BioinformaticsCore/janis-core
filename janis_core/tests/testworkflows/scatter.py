
from typing import Optional

from janis_core import (
    Workflow,
    ScatterDescription,
    ScatterMethods,
    CommandTool,
    ToolInput,
    ToolOutput,
    WildcardSelector
)
from janis_core.types import (
    String,
    File,
    Array,
    Int,
)
from janis_core.tests.testtools import (
    FileTestTool,
    StringTestTool,
    IntTestTool,
    ComponentsMandatoryTestTool,
    ArrayFileTestTool,
)

from janis_bioinformatics.data_types.bam import BamBai


# --------- #
#  SCATTER  #
# --------- #

# BasicScatterTestWF
# Scatter with no chaining

class BasicScatterTestWF(Workflow):

    def constructor(self):
        self.input('inStrArray', Array(String))
        self.input('inFileArray', Array(File))
        self.input('inIntArray', Array(Int))

        self.step(
            "stp1", 
            FileTestTool(inp=self.inFileArray),
            scatter="inp"
        )
        self.step(
            "stp2", 
            FileTestTool(inp=self.stp1.out),
            scatter="inp"
        )
        self.step(
            "stp3", 
            StringTestTool(inp=self.inStrArray),
            scatter="inp"
        )
        self.step(
            "stp4", 
            IntTestTool(inp=self.inIntArray),
            scatter="inp"
        )

        self.output("outFile", Array(File), source=self.stp1.out)
        self.output("outFileConn", Array(File), source=self.stp2.out)
        self.output("outString", Array(File), source=self.stp3.out)
        self.output("outInt", Array(File), source=self.stp4.out)

    def friendly_name(self):
        return "TEST: BasicScatterTestWF"

    def id(self) -> str:
        return self.__class__.__name__


# ChainedScatterTestWF
# Scatter with subsequent consuming steps 

class ChainedScatterTestWF(Workflow):

    def constructor(self):
        self.input('inStrArray', Array(String))
        self.input('inFileArray', Array(File))
        self.input('inIntArray', Array(Int))

        self.step(
            "stp1", 
            FileTestTool(inp=self.inFileArray),
            scatter="inp"
        )
        self.step(
            "stp2", 
            FileTestTool(inp=self.stp1.out),
            scatter="inp"
        )
        self.step(
            "stp3", 
            ArrayFileTestTool(inp=self.inFileArray),
        )
        self.step(
            "stp4", 
            FileTestTool(inp=self.stp3.out),
            scatter="inp"
        )

        self.output("outFile1", Array(File), source=self.stp1.out)
        self.output("outFile2", Array(File), source=self.stp2.out)
        self.output("outFile3", Array(File), source=self.stp3.out)
        self.output("outFile4", Array(File), source=self.stp4.out)

    def friendly_name(self):
        return "TEST: BasicScatterTestWF"

    def id(self) -> str:
        return self.__class__.__name__



# MultiFieldScatterTestWF
# Multi-field scatter (dot) with subsequent consuming steps 

class ScatterDotTestWF(Workflow):

    def constructor(self):
        self.input('inStrArray', Array(String))
        self.input('inFileArray', Array(File))
        self.input('inIntArray', Array(Int))

        self.step(
            "stp1", 
            ComponentsMandatoryTestTool(
                pos_basic=self.inFileArray,
                opt_basic=self.inStrArray
            ),
            scatter=ScatterDescription(fields=["pos_basic", "opt_basic"], method=ScatterMethods.dot)
        )
        # self.step(
        #     "stp2", 
        #     FileTestTool(
        #         inp=self.stp1.out
        #     ),
        #     scatter="inp"
        # )

        self.output("outFile1", Array(File), source=self.stp1.out)
        # self.output("outFile2", Array(File), source=self.stp2.out)

    def friendly_name(self):
        return "TEST: BasicScatterTestWF"

    def id(self) -> str:
        return self.__class__.__name__




class ScatterCrossTestWF(Workflow):

    def constructor(self):
        self.input('inStrArray', Array(String))
        self.input('inFileArray', Array(File))
        self.input('inIntArray', Array(Int))

        self.step(
            "stp1", 
            ComponentsMandatoryTestTool(
                pos_basic=self.inFileArray,
                opt_basic=self.inStrArray
            ),
            scatter=ScatterDescription(fields=["pos_basic", "opt_basic"], method=ScatterMethods.cross)
        )
        # self.step(
        #     "stp2", 
        #     FileTestTool(
        #         inp=self.stp1.out
        #     ),
        #     scatter="inp"
        # )

        self.output("outFile1", Array(File), source=self.stp1.out)
        # self.output("outFile2", Array(File), source=self.stp2.out)

    def friendly_name(self):
        return "TEST: BasicScatterTestWF"

    def id(self) -> str:
        return self.__class__.__name__




class ComprehensiveScatterTestWF(Workflow):

    def constructor(self):
        self.input('inFileArray', Array(File))
        self.input('inBamBaiArray', Array(BamBai))

        # non-secondary
        self.step(
            "prestep1", 
            FileTestTool(inp=self.inFileArray),
            scatter="inp"
        )
        self.step(
            "prestep2", 
            FileArrayTestTool(inp=self.inFileArray),
        )
        self.step(
            "scatter_to_scatter", 
            FileTestTool(inp=self.prestep1.out),
            scatter="inp"
        )
        self.step(
            "scatter_to_array", 
            FileArrayTestTool(inp=self.prestep1.out),
        )
        self.step(
            "array_to_scatter", 
            FileTestTool(inp=self.prestep2.out),
            scatter="inp"
        )

        # secondary
        self.step(
            "prestep3", 
            BamBaiTestTool(inp=self.inBamBaiArray),
            scatter="inp"
        )
        self.step(
            "scatter_secondary_to_scatter_secondary", 
            BamBaiTestTool(inp=self.prestep3.out),
            scatter="inp"
        )
        self.step(
            "scatter_secondary_to_secondary_array", 
            BamBaiArrayTestTool(inp=self.prestep3.out),
        )
        self.step(
            "secondary_array_to_scatter_secondary", 
            BamBaiTestTool(inp=self.inBamBaiArray),
            scatter="inp"
        )

    def friendly_name(self):
        return "TEST: BasicScatterTestWF"

    def id(self) -> str:
        return self.__class__.__name__


class FileArrayTestTool(CommandTool):
    def tool(self) -> str:
        return "FileArrayTestTool"

    def base_command(self) -> Optional[str | list[str]]:
        return "echo"

    def inputs(self) -> list[ToolInput]:
        return [ToolInput("inp", Array(File()))]

    def outputs(self):
        return [ToolOutput("out", Array(File()), selector=WildcardSelector('*.fasta'))]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"


class BamBaiTestTool(CommandTool):
    def tool(self) -> str:
        return "BamBaiTestTool"

    def base_command(self) -> Optional[str | list[str]]:
        return "echo"

    def inputs(self) -> list[ToolInput]:
        return [ToolInput("inp", BamBai())]

    def outputs(self):
        return [
            ToolOutput(
                "out", 
                BamBai(), 
                selector=WildcardSelector("*.bam"),
                secondaries_present_as={".bai": ".bai"},
            ),
        ]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"


class BamBaiArrayTestTool(CommandTool):
    def tool(self) -> str:
        return "BamBaiArrayTestTool"

    def base_command(self) -> Optional[str | list[str]]:
        return "echo"

    def inputs(self) -> list[ToolInput]:
        return [ToolInput("inp", Array(BamBai()))]

    def outputs(self):
        return [
            ToolOutput(
                "out", 
                BamBai(), 
                selector=WildcardSelector("*.bam"),
                secondaries_present_as={".bai": ".bai"},
            ),
        ]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"