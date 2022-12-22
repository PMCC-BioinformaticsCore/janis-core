



from janis_core import (
    Workflow,
)

from janis_core.types import (
    String,
    File,
    Array,
    Int,
    Boolean,
)
from janis_core.tests.testtools import (
    FileTestTool,
    StringTestTool,
    ComponentsTestTool,
    ArrayFileTestTool,
    ArrayComponentsTestTool,
)


# StepInputsTestWF
# calling same tool with different step inputs

class StepInputsTestWF(Workflow):

    def constructor(self):
        self.input('inFile', File)
        self.input('inStr', String)
        self.input('inInt', Int)
        self.input('inBool', Boolean)

        # full inputs
        self.step(
            "stp1", 
            ComponentsTestTool(
                pos_basic=self.inFile,
                pos_default=self.inInt,
                pos_optional=self.inStr,
                flag_true=self.inBool,
                flag_false=self.inBool,
                opt_basic=self.inStr,
                opt_default=self.inInt,
                opt_optional=self.inStr,
            )
        )
        # full inputs static
        self.step(
            "stp2", 
            ComponentsTestTool(
                pos_basic=self.inFile,
                pos_default=100,
                pos_optional="static",
                flag_true=False,
                flag_false=True,
                opt_basic="static",
                opt_default=100,
                opt_optional='',
            )
        )
        # partial inputs static
        self.step(
            "stp3", 
            ComponentsTestTool(
                pos_basic=self.inFile,
                opt_basic="static",
                opt_default=100,
                opt_optional='',
            )
        )
        # minimal inputs
        self.step(
            "stp4", 
            ComponentsTestTool(
                pos_basic=self.inFile,
                opt_basic=self.inStr,
            )
        )

        self.output("outFile1", File, source=self.stp1.out)
        self.output("outFile2", File, source=self.stp2.out)
        self.output("outFile3", File, source=self.stp3.out)
        self.output("outFile4", File, source=self.stp4.out)

    def friendly_name(self):
        return "TEST: StepInputsTestWF"

    def id(self) -> str:
        return self.__class__.__name__


class StepInputsWFInputTestWF(Workflow):

    def constructor(self):
        self.input('inFile', File)
        self.input('inFileOpt', File(optional=True))
        self.input('inStr', String)
        self.input('inInt', Int)
        self.input('inBool', Boolean)

        # full inputs
        self.step(
            "stp1", 
            ComponentsTestTool(
                pos_basic=self.inFile,
                pos_basic2=self.inFileOpt,
                pos_default=self.inInt,
                pos_optional=self.inStr,
                flag_true=self.inBool,
                flag_false=self.inBool,
                opt_basic=self.inStr,
                opt_default=self.inInt,
                opt_optional=self.inStr,
            )
        )
        self.output("outFile1", File, source=self.stp1.out)

    def friendly_name(self):
        return "TEST: StepInputsWFInputTestWF"

    def id(self) -> str:
        return self.__class__.__name__


class StepInputsStaticTestWF(Workflow):

    def constructor(self):
        self.input('inFile', File)
        self.input('inStr', String)
        self.input('inInt', Int)
        self.input('inBool', Boolean)

        # full inputs static
        self.step(
            "stp2", 
            ComponentsTestTool(
                pos_basic=self.inFile,
                pos_default=100,
                pos_optional="static",
                flag_true=False,
                flag_false=True,
                opt_basic="static",
                opt_default=100,
                opt_optional='',
            )
        )

        self.output("outFile2", File, source=self.stp2.out)

    def friendly_name(self):
        return "TEST: StepInputsStaticTestWF"

    def id(self) -> str:
        return self.__class__.__name__


class StepInputsPartialStaticTestWF(Workflow):

    def constructor(self):
        self.input('inFile', File)
        self.input('inStr', String)
        self.input('inInt', Int)
        self.input('inBool', Boolean)

        # partial inputs static
        self.step(
            "stp3", 
            ComponentsTestTool(
                pos_basic=self.inFile,
                opt_basic="static",
                opt_default=100,
                opt_optional='',
            )
        )

        self.output("outFile3", File, source=self.stp3.out)

    def friendly_name(self):
        return "TEST: StepInputsMinimalTestWF"

    def id(self) -> str:
        return self.__class__.__name__


class StepInputsMinimalTestWF(Workflow):

    def constructor(self):
        self.input('inFile', File)
        self.input('inStr', String)
        self.input('inInt', Int)
        self.input('inBool', Boolean)

        # minimal inputs
        self.step(
            "stp4", 
            ComponentsTestTool(
                pos_basic=self.inFile,
                opt_basic=self.inStr,
            )
        )
        self.output("outFile4", File, source=self.stp4.out)

    def friendly_name(self):
        return "TEST: StepInputsMinimalTestWF"

    def id(self) -> str:
        return self.__class__.__name__


# ConnectionsTestWF 
# File, String, Int workflow inputs & tool inputs 

class StepConnectionsTestWF(Workflow):

    def constructor(self):
        self.input('inStr', String)
        self.input('inFile', File)
        self.input('inInt', Int)

        self.step(
            "stp1", 
            StringTestTool(inp=self.inStr)
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



# ArrayStepInputsTestWF
# calling same tool with different step input arrays
class ArrayStepInputsTestWF(Workflow):

    def constructor(self):
        self.input('inFileArray', Array(File))
        self.input('inFileArrayOpt', Array(File, optional=True))
        self.input('inStrArray', Array(String))
        self.input('inIntArray', Array(Int))
        # self.input('inBoolArray', Array(Boolean))

        # full inputs
        self.step(
            "stp1",
            ArrayComponentsTestTool(
                pos_basic=self.inFileArray,
                pos_basic2=self.inFileArrayOpt,
                pos_default=self.inIntArray,
                pos_optional=self.inStrArray,
                # flag_true=self.inBoolArray,
                # flag_false=self.inBoolArray,
                opt_basic=self.inStrArray,
                opt_default=self.inIntArray,
                opt_optional=self.inStrArray,
            )
        )
        # full inputs static
        self.step(
            "stp2", 
            ArrayComponentsTestTool(
                pos_basic=self.inFileArray,
                pos_default=[4,5,6],
                pos_optional=["hi", "there", "friend"],
                # flag_true=[True],
                # flag_false=[True],
                opt_basic=["hi", "there", "friend"],
                opt_default=[4,5,6],
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
                inp=self.inFileArray,
            ),
        )
        self.step(
            "stp2",
            ArrayFileTestTool(
                inp=self.stp1.out,
            ),
        )
        self.output("outFiles1", source=self.stp1.out)
        self.output("outFiles2", source=self.stp2.out)

    def friendly_name(self):
        return "TEST: ArrayStepConnectionsTestWF"

    def id(self) -> str:
        return self.__class__.__name__