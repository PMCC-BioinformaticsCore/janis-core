

from janis_core import (
    Workflow,
)

from janis_core.types import (
    String,
    File,
    Array,
    Int,
    Float
)
from janis_core.tests.testtools import (
    FileTestTool,
    StringTestTool,
    IntTestTool,
    WildcardSelectorTestTool,
    FileInputSelectorTestTool,
    StringInputSelectorTestTool,
    ResourcesTestTool,
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


# directives
class DirectivesTestWF(Workflow):

    def constructor(self):
        self.input('inFile', File)

        self.step(
            "stp1", 
            ResourcesTestTool(
                inp=self.inFile,
                threads=4
            )
        )
        self.step(
            "stp2", 
            FileTestTool(inp=self.inFile)
        )

        self.output("outFile", File, source=self.stp1.out)
        self.output("outFile2", File, source=self.stp2.out)

    def friendly_name(self):
        return "TEST: DirectivesTestWF"

    def id(self) -> str:
        return self.__class__.__name__




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







