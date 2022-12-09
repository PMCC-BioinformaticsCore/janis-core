

from janis_core import (
    Workflow,
)
from janis_core.types import (
    String,
    File,
    Int,
)
from janis_core.tests.testtools import (
    FileTestTool,
    StringTestTool,
    IntTestTool,
)

class ApplesWorkflow(Workflow):
    def constructor(self):
        self.input('inStr', String)
        self.input('inInt', Int)

        self.step(
            "string_tool", 
            StringTestTool(inp=self.inStr)
        )
        self.step(
            "oranges_subworkflow", 
            OrangesWorkflow(
                inFile=self.string_tool.out,
                inInt=self.inInt
            )
        )

        self.output("outStringFile", File, source=self.string_tool.out)
        self.output("outIntFile", File, source=self.oranges_subworkflow.out)

    def friendly_name(self):
        return "TEST: ApplesWorkflow"

    def id(self) -> str:
        return self.__class__.__name__
    

class OrangesWorkflow(Workflow):
    def constructor(self):
        self.input('inFile', File)
        self.input('inInt', Int)

        self.step(
            "file_tool", 
            FileTestTool(inp=self.inFile)
        )
        self.step(
            "int_tool", 
            IntTestTool(inp=self.inInt)
        )

        self.output("out", File, source=self.int_tool.out)

    def friendly_name(self):
        return "TEST: OrangesWorkflow"

    def id(self) -> str:
        return self.__class__.__name__


class SubworkflowTestWF(Workflow):
    def constructor(self):
        self.input('inFile', File)
        self.input('inStr', String)
        self.input('inInt', Int)

        self.step(
            "file_tool", 
            FileTestTool(inp=self.inFile)
        )
        self.step(
            "apples_subworkflow", 
            ApplesWorkflow(
                inStr=self.inStr,
                inInt=self.inInt,
            )
        )

        self.output("outFile", File, source=self.file_tool.out)
        self.output("outString", File, source=self.apples_subworkflow.outStringFile)
        self.output("outInt", File, source=self.apples_subworkflow.outIntFile)

    def friendly_name(self):
        return "TEST: SubworkflowTestWF"

    def id(self) -> str:
        return self.__class__.__name__
