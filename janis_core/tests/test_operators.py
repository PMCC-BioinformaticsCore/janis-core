import unittest
from typing import List, Optional, Union

from janis_core.operators.selectors import InputSelector
from janis_core.operators.standard import BasenameOperator
from janis_core.types import Stdout, File
from janis_core.tests.testtools import SingleTestTool, EchoTestTool

from janis_core import ToolOutput, ToolInput, WorkflowBuilder
from janis_core.tool.commandtool import CommandToolBuilder, CommandTool, ToolArgument

from janis_core.operators import *
from janis_core.operators.logical import *


class TestOperators(unittest.TestCase):
    def test_not_operator(self):
        op = NotOperator(1)
        self.assertEqual("!(1)", str(op))


class TestAndOperator(unittest.TestCase):
    def test_add_operator(self):
        op = AndOperator("cond1", "cond2")
        self.assertEqual("(cond1 and cond2)", str(op))

    def test_nested_add_operator(self):
        op = AndOperator("cond1", AndOperator("cond2", "cond3"))
        self.assertEqual("(cond1 and (cond2 and cond3))", str(op))

    def test_and_two_operator(self):
        op = AndOperator("cond1", "cond2").op_and("cond3")
        self.assertEqual("((cond1 and cond2) and cond3)", str(op))


class TestAddOperator(unittest.TestCase):
    def test_add_operator(self):
        op = AddOperator(1, 2)
        self.assertEqual("(1 + 2)", str(op))

    def test_nested_add_operator(self):
        op = AddOperator(1, AddOperator(2, 3))
        self.assertEqual("(1 + (2 + 3))", str(op))

    def test_radd_to_number(self):
        op = 1 + AddOperator(2, 3)
        self.assertEqual("(1 + (2 + 3))", str(op))

    def test_add_to_number(self):
        op = AddOperator(1, 2) + 3
        self.assertEqual("((1 + 2) + 3)", str(op))


class TestBasenameOperator(unittest.TestCase):
    class TestBasenameTool(CommandTool):
        def tool(self) -> str:
            return "test_basename"

        def base_command(self) -> Optional[Union[str, List[str]]]:
            return "echo"

        def arguments(self):
            return [ToolArgument(BasenameOperator(InputSelector("inp")), position=1)]

        def inputs(self) -> List[ToolInput]:
            return [ToolInput("inp", File)]

        def outputs(self) -> List[ToolOutput]:
            return [ToolOutput("out", Stdout)]

        def container(self) -> str:
            return "ubuntu:latest"

        def version(self) -> str:
            return None

    def test_cwl(self):
        TestBasenameOperator.TestBasenameTool().translate("wdl")


class TestIfOrElseOperator(unittest.TestCase):
    def test_if(self):
        wf = WorkflowBuilder("wf")
        wf.input("inp", int)
        op = If(wf.inp > 5, wf.inp, -1)
        self.assertEqual("(inputs.inp > 5) ? inputs.inp : -1", str(op))


class TestIndexOperator(unittest.TestCase):
    @classmethod
    def setUpClass(cls):

        w = WorkflowBuilder("test_operators")

        w.input("inp", Array(File()))
        inval = w.inp[0].basename()
        w.step("echo", SingleTestTool(input1=inval))
        w.output("out", source=w.echo)
        cls.wf = w

        w2 = WorkflowBuilder("test_scattered_operator_with_alias")

        w2.input("inp", Array(Array(String)))
        w2.step("echo", SingleTestTool(input1=w2.inp[0]), scatter="input1")
        w2.output("out", source=w.echo)
        cls.wf2 = w2

    def test_wdl(self):
        self.wf.translate("wdl", allow_empty_container=True)

    # def test2_wdl(self):
    #     self.assertRaises(
    #         Exception, self.wf2.translate, translation="wdl", allow_empty_container=True
    #     )


class TestMixedOperators(unittest.TestCase):
    def test_expression_default(self):

        wf = WorkflowBuilder("test_expression_defaults")
        wf.input("inp", Optional[str])

        wf.step(
            "echo",
            EchoTestTool(inp="Hello, " + If(IsDefined(wf.inp), wf.inp, ", Michael!")),
        )

        wf.output("out", source=wf.echo)

        wf.translate("cwl")
