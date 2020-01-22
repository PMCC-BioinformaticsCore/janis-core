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

        w.step(
            "cat",
            Cat(file=w.name),
            when=AndOperator(w.echo.out == "Hello, Michael", True),
        )

        w.output("out", source=w.echo.out)

        # w.translate(
        #     "wdl", to_disk=True, export_path="~/Desktop/tmp/{name}", validate=True
        # )
