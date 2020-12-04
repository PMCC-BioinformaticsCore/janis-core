import os
import operator
from typing import Optional, List, Dict, Union
from unittest import TestCase
from janis_core.tool.test_suite_runner import ToolTestSuiteRunner
from janis_core.tool.test_classes import (
    TTestCase,
    TTestExpectedOutput,
    TTestPreprocessor,
)


from janis_core import ToolOutput, ToolInput, CommandTool, Stdout, File


class TestTool(CommandTool):
    def tool(self) -> str:
        return "EchoTestTool"

    def base_command(self) -> Optional[Union[str, List[str]]]:
        return "echo"

    def inputs(self) -> List[ToolInput]:
        return [ToolInput("inp", str, position=0)]

    def outputs(self):
        return [ToolOutput("out", Stdout), ToolOutput("file", Stdout)]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"


class TestToolTestRunner(TestCase):
    def setUp(self):
        self.tool = TestTool()
        self.test_data_dir = os.path.join(os.getcwd(), "data")

    def test_get_expected_value(self):
        runner = ToolTestSuiteRunner(self.tool)
        t1 = TTestExpectedOutput(
            tag="out",
            preprocessor=TTestPreprocessor.Value,
            operator=operator.eq,
            expected_value=5,
        )
        assert runner.get_expected_value(t1) == 5

        t2 = TTestExpectedOutput(
            tag="out",
            preprocessor=TTestPreprocessor.Value,
            operator=operator.eq,
            expected_value=["abc", "def"],
        )
        assert runner.get_expected_value(t2) == ["abc", "def"]

    def test_get_value_to_compare(self):
        runner = ToolTestSuiteRunner(self.tool)

        t1 = TTestExpectedOutput(
            tag="out",
            preprocessor=TTestPreprocessor.Value,
            operator=operator.eq,
            expected_value="xxx",
        )

        assert runner.get_value_to_compare(t1, 5) == 5

        t2 = TTestExpectedOutput(
            tag="file",
            preprocessor=TTestPreprocessor.FileContent,
            operator=operator.eq,
            expected_value="xxx",
        )

        file_path = os.path.join(self.test_data_dir, "test.txt")
        assert runner.get_value_to_compare(t2, file_path) == "abc\ndef\nghi\n"

        t2 = TTestExpectedOutput(
            tag="file",
            preprocessor=TTestPreprocessor.LineCount,
            operator=operator.eq,
            expected_value="xxx",
        )

        file_path = os.path.join(self.test_data_dir, "test.txt")
        assert runner.get_value_to_compare(t2, file_path) == 3
