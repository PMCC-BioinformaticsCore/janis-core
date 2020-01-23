import unittest

from janis_core.types import InputSelector, AndOperator, Operator
from janis_unix import Echo, Cat

from janis_core.workflow.workflow import WorkflowBuilder


class TestConditionals(unittest.TestCase):
    def test_1(self):
        w = WorkflowBuilder("conditionalTest")

        w.input("inp", int, value=1)
        w.input("name", str, value="Michael")

        cond = w.name.as_operator() == "Michael"
        w.step("echo", Echo(inp=w.name), when=cond)

        w.step("cat", Cat(file=w.echo.out), when=w.echo.out == "Hello, Michael")

        w.output("out", source=w.echo.out)

        w.translate(
            "wdl", to_disk=True, export_path="~/Desktop/tmp/{name}", validate=True
        )

    def test_switch(self):

        w = WorkflowBuilder("switchTest")

        w.input("inp", int, value=2)
        w.input("inp1", str, value="Hello")
        w.input("inp2", str, value="Hi there")

        w.switch(
            "echoswitch",
            [(w.inp.as_operator() > 1, Echo(inp=w.inp1)), Echo(inp=w.inp2)],
        )

        w.output("out", source=w.echoswitch)

        w.translate(
            "wdl", to_disk=True, export_path="~/Desktop/tmp/{name}", validate=True
        )
