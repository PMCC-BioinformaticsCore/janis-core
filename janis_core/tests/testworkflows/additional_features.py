

from typing import Optional
from janis_core.operators.logical import If, IsDefined
from janis_core import (
    String,
    Array,
    File,
    Workflow,
    FirstOperator,
    ForEachSelector,
)

from janis_core.tests.testtools import (
    EchoTestTool,
    CatTestTool,
    SecondaryOutputTestTool,
)



class StepInputExpressionTestWF(Workflow):
    def constructor(self):
        self.input("mystring", Optional[str], value="")
        self.input("mystring_backup", Optional[str])

        self.step(
            "print",
            EchoTestTool(
                inp=If(IsDefined(self.mystring), self.mystring, self.mystring_backup)
            ),
        )

        self.output("out", source=self.print)

    def friendly_name(self):
        return "TEST: WorkflowWithStepInputExpression"

    def id(self) -> str:
        return self.__class__.__name__
        
        
class ArraysOfSecondaryFilesOutputsTestWF(Workflow):
    def id(self) -> str:
        return "ArraysOfSecondaryFilesOutputsTestWF"

    def friendly_name(self):
        return "Test Workflow That outputs arrays of secondary files"

    def constructor(self):
        self.input("inp", Array(String))

        self.step(
            "stp", SecondaryOutputTestTool(testtool=self.inp), scatter="testtool"
        )

        self.output("out", source=self.stp.out)


class ConditionStepTestWF(Workflow):
    def constructor(self):
        self.input("mystring", Optional[str], value=None)

        someString = FirstOperator(
            [
                self.mystring,
                self.step(
                    "get_string",
                    EchoTestTool(inp="Some default value"),
                    when=self.mystring.is_null(),
                ).out,
            ]
        )

        self.step(
            "print",
            EchoTestTool(inp=someString),
        )

        self.output("out", source=self.print)

    def friendly_name(self):
        return "TEST: TestWorkflowWithConditionStep"

    def id(self) -> str:
        return self.__class__.__name__


class AliasSelectorTestWF(Workflow):
    def id(self) -> str:
        return "TestWorkflowWithAliasSelectorWorkflow"

    def friendly_name(self):
        return "Test Workflow with alias selector in the output"

    def constructor(self):
        self.input("inp", String, value="abc")

        self.step("stp1", SecondaryOutputTestTool(testtool=self.inp))
        self.step("stp2", CatTestTool(inp=self.stp1.out.as_type(File)))

        self.output("out", source=self.stp1.out)


class ForEachTestWF(Workflow):
    def constructor(self):
        self.input("inp", Array(str))
        self.step(
            "print", EchoTestTool(inp=ForEachSelector() + "-hello"), _foreach=self.inp
        )
        self.output("out", source=self.print.out)

    def friendly_name(self):
        return self.id()

    def id(self) -> str:
        return "ForEachTestWF"

