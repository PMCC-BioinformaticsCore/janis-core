
from typing import Optional
from janis_core import (
    Workflow,
    ToolInput,
    ToolArgument,
    ToolOutput,
    CommandToolBuilder,
    WorkflowBuilder,
    CeilOperator,
    FileSizeOperator,
    MultiplyOperator,
    DivideOperator,
    InputSelector,
    AddOperator,
)



from janis_core.types import (
    String,
    File,
    Array,
    Int,
    Float,
    Stdout
)
from janis_core.tests.testtools import (
    FileTestTool,
    StringTestTool,
    IntTestTool,
    WildcardSelectorTestTool,
    FileInputSelectorTestTool,
    StringInputSelectorTestTool,
    ArrayFileTestTool,
    ArrayIntTestTool,
    ArrayStringTestTool,
    ArrayWildcardSelectorTestTool,
    ArrayInputSelectorTestTool,
)

# ------------- #
#  BASIC TYPES  #
# ------------- #

# BasicInTypesTestWF 
# File, String, Int workflow inputs & tool inputs 

class BasicIOTestWF(Workflow):

    def constructor(self):
        self.input('inFile', File)
        self.input('inStr', String)
        self.input('inInt', Int)

        self.step(
            "stp1", 
            FileTestTool(inp=self.inFile)
        )
        self.step(
            "stp2", 
            StringTestTool(inp=self.inStr)
        )
        self.step(
            "stp3", 
            IntTestTool(inp=50)
        )
        # self.step(
        #     "stp4", 
        #     FileTestTool(inp=self.stp1.out)
        # )

        self.output("outFile", File, source=self.stp1.out)
        self.output("outString", File, source=self.stp2.out)
        self.output("outInt", File, source=self.stp3.out)
        # self.output("outFileConn", File, source=self.stp4.out)

    def friendly_name(self):
        return "TEST: BasicIOTestWF"

    def id(self) -> str:
        return self.__class__.__name__


# all WildcardSelector use cases
class WildcardSelectorOutputTestWF(Workflow):

    def constructor(self):
        self.input('inFile', File)
        self.input('inFileArr', Array(File))

        self.step(
            "stp1", 
            WildcardSelectorTestTool(inp=self.inFile)
        )
        self.step(
            "stp2", 
            ArrayWildcardSelectorTestTool(inp=self.inFileArr)
        )

        self.output("outFile", File, source=self.stp1.out)
        self.output("outFileArr", Array(File), source=self.stp2.out)

    def friendly_name(self):
        return "TEST: WildcardSelectorOutputTestWF"

    def id(self) -> str:
        return self.__class__.__name__


# all InputSelector use cases
class InputSelectorTestWF(Workflow):

    def constructor(self):
        self.input('inFile', File)
        self.input('inStr', String)
        self.input('inFileArr', Array(File))

        self.step(
            "stp1", 
            FileInputSelectorTestTool(inp=self.inFile)
        )
        self.step(
            "stp2", 
            StringInputSelectorTestTool(inp=self.inStr)
        )
        self.step(
            "stp3", 
            ArrayInputSelectorTestTool(inp=self.inFileArr)
        )

        self.output("outFile", File, source=self.stp1.out)
        self.output("outFileArr", File, source=self.stp2.out)

    def friendly_name(self):
        return "TEST: InputSelectorOutputTestWF"

    def id(self) -> str:
        return self.__class__.__name__


# DIRECTIVES
resources_tool = CommandToolBuilder(
    tool="ResourcesTestToolCB",
    version="TEST",
    base_command="echo",
    inputs=[
        ToolInput("inp", File, position=1), 
        ToolInput("threads", Int)
    ],
    outputs=[ToolOutput("out", Stdout)],
    container="quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0",
    disk=CeilOperator(
            AddOperator(
                DivideOperator(
                    FileSizeOperator(InputSelector('inp')), 
                    MultiplyOperator(MultiplyOperator(1024, 1024), 1024)), 
                20)),
    memory=MultiplyOperator(15, 1024),
    cpus=1,
    time=60,
)

direcs_wf = WorkflowBuilder("DirectivesTestWF")
direcs_wf.input("inFile", File)
direcs_wf.step("stp1", resources_tool(inp=direcs_wf.inFile, threads=4))
direcs_wf.output("outFile", source=direcs_wf.stp1.out)
DirectivesTestWF = direcs_wf
    

# -------- #
#  ARRAYS  #
# -------- #

# BasicArrayInOutTypesTestWF 
# File, String, Int workflow & tool inputs / outputs

class ArrayIOTestWF(Workflow):

    def constructor(self):
        self.input('inFileArray', Array(File))
        self.input('inStrArray', Array(String))
        self.input('inIntArray', Array(Int))
        self.step(
            "stp1",
            ArrayStringTestTool(
                inp=self.inStrArray,
            ),
        )
        self.step(
            "stp2",
            ArrayFileTestTool(
                inp=self.inFileArray,
            ),
        )
        self.step(
            "stp3",
            ArrayIntTestTool(
                inp=self.inIntArray,
            ),
        )
        self.output("outStrings", source=self.stp1.out)
        self.output("outFiles", source=self.stp2.out)
        self.output("outInts", source=self.stp3.out)

    def friendly_name(self):
        return "TEST: ArrayIOTestWF"

    def id(self) -> str:
        return self.__class__.__name__



class ArrayIOExtrasTestWF(Workflow):

    def constructor(self):
        self.input('inFileArray', Array(File))
        self.input('inStrArray', Array(String))
        self.input('inIntArray', Array(Int))
        self.input('inFloatArray', Array(Float))
        self.step(
            "stp1",
            ArrayStringTestTool(
                inp=self.inStrArray,
            ),
        )
        self.step(
            "stp2",
            ArrayFileTestTool(
                inp=self.inFileArray,
            ),
        )
        self.step(
            "stp3",
            ArrayIntTestTool(
                inp=self.inIntArray,
            ),
        )
        self.step(
            "stp4",
            ArrayStringTestTool(
                inp=['hello', 'there!'],
            ),
        )
        self.output("outStrings", source=self.stp1.out)
        self.output("outFiles", source=self.stp2.out)
        self.output("outInts", source=self.stp3.out)
        self.output("outStrings2", source=self.stp4.out)

    def friendly_name(self):
        return "TEST: ArrayIOTestWF"

    def id(self) -> str:
        return self.__class__.__name__







