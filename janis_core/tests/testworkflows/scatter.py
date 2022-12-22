
from janis_core import (
    Workflow,
    ScatterDescription,
    ScatterMethods
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
    ComponentsTestTool,
    ArrayFileTestTool,
)


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
            ComponentsTestTool(
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
            ComponentsTestTool(
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



