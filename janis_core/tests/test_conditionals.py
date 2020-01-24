import unittest

from janis_core.utils import is_module_available
from janis_core.workflow.workflow import WorkflowBuilder


@unittest.skipUnless(is_module_available("janis_unix"), "janis_unix is not available")
class TestConditionals(unittest.TestCase):
    def test_1(self):
        from janis_unix import Echo, Cat

        w = WorkflowBuilder("conditionalTest")

        w.input("inp", int, value=1)
        w.input("name", str, value="Michael")

        w.step("echo", Echo(inp=w.name), when=w.inp > 1)
        w.step("cat", Cat(file=w.echo.out), when=w.echo.out == "Hello, Michael")

        w.output("out", source=w.echo.out)

        w.translate(
            "wdl"  # to_disk=True, export_path="~/Desktop/tmp/{name}", validate=True
        )

    def test_switch(self):
        from janis_unix import Echo, Cat

        w = WorkflowBuilder("switchTest")

        w.input("inp", int, value=2)
        w.input("inp1", str, value="Hello")
        w.input("inp2", str, value="Hi there")

        w.conditional("echoswitch", [(w.inp > 1, Echo(inp=w.inp1)), Echo(inp=w.inp2)])

        w.output("out", source=w.echoswitch)

        w.translate(
            "wdl"  # to_disk=True, export_path="~/Desktop/tmp/{name}", validate=True
        )
