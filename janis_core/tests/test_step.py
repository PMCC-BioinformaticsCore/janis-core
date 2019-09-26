from unittest import TestCase

from janis_core import String, CommandTool, ToolInput, Array, Logger, WorkflowBuilder
from janis_core.graph.steptaginput import StepTagInput, Edge
from janis_core.tests.testtools import SingleTestTool


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
        self.wf = WorkflowBuilder("testStep")

    def tearDown(self):
        Logger.unmute()

    def test_sources_single_slashed(self):
        start1 = self.wf.input("test1", str)
        step = self.wf.step("step", SingleTestTool(inputs=start1))
        source = step.sources["inputs"].slashed_source()

        self.assertEqual(start1.id(), source)

    def test_sources_single_dotted(self):
        start1 = self.wf.input("test1", str)
        step = self.wf.step("step", SingleTestTool(inputs=start1))
        source = step.sources["inputs"].dotted_source()

        self.assertEqual(start1.id(), source)

    def test_sources_multiple_slashed(self):
        test1 = self.wf.input("test1", str)
        test2 = self.wf.input("test2", str)

        step = self.wf.step("step", ArrayTestTool(inputs=[test1, test2]))
        source = step.sources["inputs"].slashed_source()

        self.assertIsInstance(source, list)
        self.assertEqual(len(source), 2)

    def test_sources_multiple_dotted(self):
        test1 = self.wf.input("test1", str)
        test2 = self.wf.input("test2", str)

        step = self.wf.step("step", ArrayTestTool(inputs=[test1, test2]))
        source = step.sources["inputs"].dotted_source()

        self.assertIsInstance(source, list)
        self.assertEqual(len(source), 2)
