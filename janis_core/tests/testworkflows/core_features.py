

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
    Boolean,
)
from janis_bioinformatics.data_types.bam import BamBai
from janis_core.tests.testtools import (
    FileTestTool,
    StringTestTool,
    IntTestTool,
    ComponentsTestTool,
    SecondariesTestTool,

    ArrayFileTestTool,
    ArrayIntTestTool,
    ArrayStringTestTool,

    ArrayComponentsTestTool,
    ArraySecondariesTestTool,
)



# ------------- #
#  BASIC TYPES  #
# ------------- #

# BasicInTypesTestWF 
# File, String, Int workflow inputs & tool inputs 

class BasicIOTestWF(Workflow):

    def constructor(self):
        self.input('inFile', File)
        self.input('inString', String)
        self.input('inInt', Int)

        self.step(
            "stp1", 
            FileTestTool(inp=self.inFile)
        )
        self.step(
            "stp2", 
            StringTestTool(inp=self.inString)
        )
        self.step(
            "stp3", 
            IntTestTool(inp=self.inInt)
        )

        self.output("outFile", File, source=self.stp1.out)
        self.output("outString", File, source=self.stp1.out)
        self.output("outInt", File, source=self.stp3.out)

    def friendly_name(self):
        return "TEST: BasicIOTestWF"

    def id(self) -> str:
        return self.__class__.__name__


# StepInputsTestWF
# calling same tool with different step inputs

class StepInputsTestWF(Workflow):

    def constructor(self):
        self.input('inFile', File)
        self.input('inString', String)
        self.input('inInt', Int)
        self.input('inBool', Boolean)

        # full inputs
        self.step(
            "stp1", 
            ComponentsTestTool(
                pos_basic=self.inFile,
                pos_default=self.inString,
                pos_optional=self.inString,
                flag_true=self.inBool,
                flag_false=self.inBool,
                opt_basic=self.inString,
                opt_default=self.inString,
                opt_optional=self.inString,
            )
        )
        # full inputs static
        self.step(
            "stp2", 
            ComponentsTestTool(
                pos_basic=self.inFile,
                pos_default="static",
                pos_optional="static",
                flag_true=False,
                flag_false=True,
                opt_basic="static",
                opt_default="static",
                opt_optional="static",
            )
        )
        # minimal inputs
        self.step(
            "stp3", 
            ComponentsTestTool(
                pos_basic=self.inFile,
                opt_basic=self.inString,
            )
        )

        self.output("outFile1", File, source=self.stp1.out)
        self.output("outFile2", File, source=self.stp2.out)
        self.output("outFile3", File, source=self.stp3.out)

    def friendly_name(self):
        return "TEST: StepInputsTestWF"

    def id(self) -> str:
        return self.__class__.__name__


# ConnectionsTestWF 
# File, String, Int workflow inputs & tool inputs 

class StepConnectionsTestWF(Workflow):

    def constructor(self):
        self.input('inString', String)
        self.input('inFile', File)
        self.input('inInt', Int)

        self.step(
            "stp1", 
            StringTestTool(inp=self.inString)
        )
        self.step(
            "stp2", 
            FileTestTool(inp=self.stp1.out)
        )

        self.output("outString", File, source=self.stp1.out)
        self.output("outFile", File, source=self.stp2.out)

    def friendly_name(self):
        return "TEST: StepConnectionsTestWF"

    def id(self) -> str:
        return self.__class__.__name__



# -------- #
#  ARRAYS  #
# -------- #

# BasicArrayInOutTypesTestWF 
# File, String, Int workflow & tool inputs / outputs

class ArrayIOTestWF(Workflow):

    def constructor(self):
        self.input('inStrArray', Array(String))
        self.input('inFileArray', Array(File))
        self.input('inIntArray', Array(Int))
        self.step(
            "stp1",
            ArrayStringTestTool(
                ins=self.inStrArray,
            ),
        )
        self.step(
            "stp2",
            ArrayFileTestTool(
                ins=self.inFileArray,
            ),
        )
        self.step(
            "stp3",
            ArrayIntTestTool(
                ins=self.inIntArray,
            ),
        )
        self.output("outStrings", source=self.stp1.outs)
        self.output("outFiles", source=self.stp2.outs)
        self.output("outInts", source=self.stp3.outs)

    def friendly_name(self):
        return "TEST: ArrayIOTestWF"

    def id(self) -> str:
        return self.__class__.__name__



# ArrayStepInputsTestWF
# calling same tool with different step input arrays
class ArrayStepInputsTestWF(Workflow):

    def constructor(self):
        self.input('inFileArray', Array(File))
        self.input('inStrArray', Array(String))
        self.input('inIntArray', Array(Int))
        self.input('inBoolArray', Array(Boolean))

        # full inputs
        self.step(
            "stp1",
            ArrayComponentsTestTool(
                pos_basic=self.inFileArray,
                pos_default=self.inStrArray,
                pos_optional=self.inStrArray,
                flag_true=self.inBoolArray,
                flag_false=self.inBoolArray,
                opt_basic=self.inStrArray,
                opt_default=self.inStrArray,
                opt_optional=self.inStrArray,
            )
        )
        # full inputs static
        self.step(
            "stp2", 
            ArrayComponentsTestTool(
                pos_basic=self.inFileArray,
                pos_default=["hi", "there", "friend"],
                pos_optional=["hi", "there", "friend"],
                flag_true=[True],
                flag_false=[True],
                opt_basic=["hi", "there", "friend"],
                opt_default=["hi", "there", "friend"],
                opt_optional=["hi", "there", "friend"],
            )
        )
        # minimal inputs
        self.step(
            "stp3", 
            ArrayComponentsTestTool(
                pos_basic=self.inFileArray,
                opt_basic=self.inStrArray,
            )
        )

        self.output("outFile1", File, source=self.stp1.out)
        self.output("outFile2", File, source=self.stp2.out)
        self.output("outFile3", File, source=self.stp3.out)

    def friendly_name(self):
        return "TEST: ArrayStepInputsTestWF"

    def id(self) -> str:
        return self.__class__.__name__


# ArrayStepConnectionsTestWF 
# Array(File) step io connections
# janis does not allow non-file outputs, so can't test
# Array(String) -> Array(String) step connections. 

class ArrayStepConnectionsTestWF(Workflow):

    def constructor(self):
        self.input('inStrArray', Array(String))
        self.input('inFileArray', Array(File))
        self.input('inIntArray', Array(Int))

        self.step(
            "stp1",
            ArrayFileTestTool(
                ins=self.inFileArray,
            ),
        )
        self.step(
            "stp2",
            ArrayFileTestTool(
                ins=self.stp1.outs,
            ),
        )
        self.output("outFiles1", source=self.stp1.outs)
        self.output("outFiles2", source=self.stp2.outs)

    def friendly_name(self):
        return "TEST: ArrayStepConnectionsTestWF"

    def id(self) -> str:
        return self.__class__.__name__





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
            StringTestTool(inp=self.inStrArray),
            scatter="inp"
        )
        self.step(
            "stp3", 
            IntTestTool(inp=self.inIntArray),
            scatter="inp"
        )

        self.output("outFile", Array(File), source=self.stp1.out)
        self.output("outString", Array(File), source=self.stp2.out)
        self.output("outInt", Array(File), source=self.stp3.out)

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

        self.output("outFile1", Array(File), source=self.stp1.out)
        self.output("outFile2", Array(File), source=self.stp2.out)

    def friendly_name(self):
        return "TEST: BasicScatterTestWF"

    def id(self) -> str:
        return self.__class__.__name__


# MultiFieldScatterTestWF
# Multi-field scatter (dot) with subsequent consuming steps 

class MultiFieldScatterTestWF(Workflow):

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
        self.step(
            "stp2", 
            FileTestTool(
                inp=self.stp1.out
            ),
            scatter="inp"
        )

        self.output("outFile1", Array(File), source=self.stp1.out)
        self.output("outFile2", Array(File), source=self.stp2.out)

    def friendly_name(self):
        return "TEST: BasicScatterTestWF"

    def id(self) -> str:
        return self.__class__.__name__



# ------------- #
#  SECONDARIES  #
# ------------- #


class SecondariesIOTestWF(Workflow):
    def id(self) -> str:
        return "SecondaryFileScatterTestWF"

    def friendly_name(self):
        return "WF which uses SecondaryFile types for workflow / tool IO"

    def constructor(self):
        self.input('inBamBai', BamBai)

        self.step(
            "stp1", 
            SecondariesTestTool(
                inp=self.inBamBai
            ), 
        )

        self.output("outBamBai", source=self.stp1.out)


# ------------------- #
#  DISGUSTING COMBOS  #
# ------------------- #


# class ArrayScatterTestWF(Workflow):
#     def id(self) -> str:
#         return "ArrayScatterTestWF"

#     def friendly_name(self):
#         return "WF which uses Array(File) and Scatter"

#     def constructor(self):
#         self.input('inStrArray', Array(String))
#         self.input('inFileArray', Array(File))
#         self.input('inIntArray', Array(Int))

#         self.step(
#             "stp1", 
#             (testtool=self.inStrArray), 
#             scatter="testtool"
#         )

#         self.output("outStrArray", source=self.stp1.out)


# TODO FAILS secondaries_present_as
class ArraySecondariesTestWF(Workflow):
    def id(self) -> str:
        return "ArraySecondariesTestWF"

    def friendly_name(self):
        return "WF which uses Arrays pf SecondaryFile types for workflow / tool IO"

    def constructor(self):
        self.input('inBamBaiArray', Array(BamBai))
        
        self.step(
            "stp1", 
            ArraySecondariesTestTool(
                inp=self.inBamBaiArray
            ), 
        )

        self.output("outBamBaiArray", source=self.stp1.outArray)
        self.output("outStdout", source=self.stp1.outStdout)


# class ScatterSecondaries(Workflow):
#     def id(self) -> str:
#         return "ScatterSecondaries"

#     def friendly_name(self):
#         return "WF which uses Scatter and Secondaries"

#     def constructor(self):
#         self.input('inStrArray', Array(String))
#         self.input('inFileArray', Array(File))
#         self.input('inIntArray', Array(Int))

#         self.step(
#             "stp1", 
#             (testtool=self.inStrArray), 
#             scatter="testtool"
#         )

#         self.output("outStrArray", source=self.stp1.out)


# class HolyGrail(Workflow):
#     def id(self) -> str:
#         return "HolyGrail"

#     def friendly_name(self):
#         return "WF which uses almost every core feature"

#     def constructor(self):
#         pass


