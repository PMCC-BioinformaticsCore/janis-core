
## type: ignore

from janis_core import Workflow
from janis_core.tests.testtools import FileTypeTestTool
from janis_core.tests.testtools import StringTypeTestTool
from janis_core.tests.testtools import IntTypeTestTool

from janis_core.types import (
    File,
    String,
    Int
)

class CoreTypesTestWF(Workflow):
    def constructor(self):
        self.input('inFile', File)
        self.input('inString', String)
        self.input('inInt', Int)

        self.step(
            "stp1", 
            FileTypeTestTool(inp=self.inFile)
        )
        self.step(
            "stp2", 
            StringTypeTestTool(inp=self.inString)
        )
        self.step(
            "stp3", 
            IntTypeTestTool(inp=self.inInt)
        )

        self.output("outFile", File, source=self.stp1.out)
        self.output("outString", File, source=self.stp1.out)
        self.output("outInt", File, source=self.stp3.out)

    def friendly_name(self):
        return "TEST: CoreTypesTestWF"

    def id(self) -> str:
        return self.__class__.__name__


wf = CoreTypesTestWF()

