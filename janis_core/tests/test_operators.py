import unittest
from typing import List, Optional, Union

from janis_core.operators.selectors import InputSelector
from janis_core.operators.standard import BasenameOperator
from janis_core.types import Stdout, File
from janis_core.tests.testtools import SingleTestTool

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
        self.assertEqual("(cond1 && cond2)", str(op))

    def test_nested_add_operator(self):
        op = AndOperator("cond1", AndOperator("cond2", "cond3"))
        self.assertEqual("(cond1 && (cond2 && cond3))", str(op))

    def test_and_to_operator(self):
        op = AndOperator("cond1", "cond2").op_and("cond3")
        self.assertEqual("((cond1 && cond2) && cond3)", str(op))


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


class TestIndexOperator(unittest.TestCase):
    @classmethod
    def setUpClass(cls):

        w = WorkflowBuilder("test_operators")

        w.input("inp", Array(String))
        w.step("echo", SingleTestTool(inputs=w.inp[0]))
        w.output("out", source=w.echo)
        cls.wf = w

    def test_wdl(self):
        self.wf.translate("wdl", allow_empty_container=True)
