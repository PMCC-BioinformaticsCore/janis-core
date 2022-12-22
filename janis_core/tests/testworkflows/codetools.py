

from janis_core import (
    Workflow,
)

from ..testtools.types import SecondaryTestType
from janis_core.types import (
    String,
    File,
    Array,
    Int,
)

from janis_core.tests.testtools import (
    FileInputPythonTestTool,
    FileOutputPythonTestTool,
    SplitTextPythonTestTool,
    JoinArrayPythonTestTool,
    SecondaryInputPythonTestTool,
    MultiTypesInputPythonTool
)


class InputsPythonToolTestWF(Workflow):

    def constructor(self):
        self.input('inFile', File)
        self.input('inSecondaryType', SecondaryTestType)
        self.input('inStr', String)
        self.input('inStrArr', Array(String))
        self.input('inInt', Int)

        self.step(
            "stp0", 
            MultiTypesInputPythonTool(
                inp1=self.inFile,
                inp2=self.inStr,
                inp3=self.inInt,
            )
        )
        self.step(
            "stp1", 
            JoinArrayPythonTestTool(inp=self.inStrArr)
        )
        self.step(
            "stp2", 
            SecondaryInputPythonTestTool(inp=self.inSecondaryType)
        )

    def friendly_name(self):
        return "TEST: PythonToolTestWF"

    def id(self) -> str:
        return self.__class__.__name__


class OutputsPythonToolTestWF(Workflow):

    def constructor(self):
        self.input('inFile', File)

        self.step(
            "stp0", 
            FileOutputPythonTestTool(inp=self.inFile)
        )
        self.step(
            "stp1", 
            FileInputPythonTestTool(inp=self.stp0.out)
        )
        self.step(
            "stp2", 
            SplitTextPythonTestTool(inp=self.stp1.out)
        )

        self.output("out0", File, source=self.stp0.out)
        self.output("out1", String, source=self.stp1.out)
        self.output("out2", Array(String), source=self.stp2.out)

    def friendly_name(self):
        return "TEST: PythonToolTestWF"

    def id(self) -> str:
        return self.__class__.__name__