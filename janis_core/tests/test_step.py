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
        return [
            ToolInput("inputs", Array(String())),
            ToolInput("optionalInput", String(optional=True)),
            ToolInput("optionalInput2", String, default="2"),
        ]

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
