from unittest import TestCase

from janis_core import String, CommandTool, ToolInput, Array, Logger, Workflow
from janis_core.graph.stepinput import StepInput, Edge


class ArrayTestTool(CommandTool):
    @staticmethod
    def tool():
        return "TestStepTool"

    @staticmethod
    def base_command():
        return "echo"

    def inputs(self):
        return [ToolInput("inputs", Array(String()))]

    def outputs(self):
        return []

    def friendly_name(self):
        return "'test_step' array of strings tool"

    @staticmethod
    def container():
        return None

    @staticmethod
    def version():
        return None


class TestStep(TestCase):
    def setUp(self):
        Logger.mute()
        self.wf = Workflow("testStep")

    def tearDown(self):
        Logger.unmute()

    def test_sources_single(self):
        start1 = self.wf.input("test1", str)
        step = self.wf.step("step", ArrayTestTool, inputs=start1)
        source = step.sources["inputs"].dotted_source()

        self.assertEqual(source, start1.id())

    def test_sources_multiple(self):
        test1 = self.wf.input("test1", str)
        test2 = self.wf.input("test2", str)

        step = self.wf.step("step", ArrayTestTool, inputs=[test1, test2])
        source = step.sources["inputs"].source()

        self.assertIsInstance(source, list)
        self.assertEqual(len(source), 2)
