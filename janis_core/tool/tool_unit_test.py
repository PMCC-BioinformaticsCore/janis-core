import unittest
import hashlib
import os
from typing import Dict, List, Optional, Any
from pkg_resources import parse_version

from janis_core.tool.tool import Tool, TTestExpectedOutput, TTestCompared, TTestCase
from janis_core.types import DataType, File, String


class ToolTestSuiteRunner():
    """
    A class to read tool test cases, run the test cases, and assert the expected output
    """
    def __init__(self, tool: Tool):
        self.tool = tool

    def run(self, input: Dict[str, str]):
        janis_assistant_version_required_min = "0.10.6"

        # Note, we do not want janis_core to depend on janis_assistant except for when we run these tests
        # So, we only import this package here in this function
        # Make sure the correct version of janis_assistant is installed first
        try:
            import janis_assistant
            if parse_version(janis_assistant.__version__) < parse_version(janis_assistant_version_required_min):
                raise Exception()

            from janis_assistant.main import run_with_outputs
        except Exception as e:
            raise Exception(f"to run this test, janis_asisstant >= {janis_assistant_version_required_min}"
                            f" must be installed")

        # Now all good, we run the test
        return run_with_outputs(self.tool, input)

    def run_all_tests(self):
        """
        Run all test cases for one tool
        """
        failed = set()
        succeeded = set()
        for t in self.tool.tests() or []:
            output = self.run(t.input)

            for expected_output in t.output:
                output_value = output[expected_output.tag]
                if self.assert_expected_output(expected_output, output_value) is False:
                    failed.add(str(expected_output))
                else:
                    succeeded.add(str(expected_output))

        return failed, succeeded

    def run_one_test_case(self, t: TTestCase):
        """
        Run one test case and assert multiple expected output
        """
        failed = set()
        succeeded = set()

        output = self.run(t.input)
        for expected_output in t.output:
            output_value = output[expected_output.tag]
            if self.assert_expected_output(expected_output, output_value) is False:
                failed.add(str(expected_output))
            else:
                succeeded.add(str(expected_output))

        return failed, succeeded

    def assert_expected_output(self, expected_output: TTestExpectedOutput, output_value: Any):
        value = self.get_value_to_compare(expected_output, output_value)
        return expected_output.operator(value, expected_output.expected_value)

    def get_value_to_compare(self, expected_output: TTestExpectedOutput, output_value: Any) -> Any:
        output_tag = expected_output.tag
        value = None

        if expected_output.compared == TTestCompared.Value:
            value = output_value
        elif expected_output.compared == TTestCompared.FileContent:
            with open(output_value, "rb") as f:
                value = f.read()
        elif expected_output.compared == TTestCompared.FileMd5:
            value = self.read_md5(output_value)
        elif expected_output.compared == TTestCompared.FileSize:
            value = os.path.getsize(output_value)
        elif expected_output.compared == TTestCompared.LineCount:
            output_type = self.tool.outputs_map().get(output_tag).outtype

            try:
                if isinstance(output_type, File):
                    # text file only here
                    value = sum(1 for line in open(output_value, "r"))
                elif isinstance(output_type, String):
                    value = len(output_value.splitlines())
            except Exception as e:
                raise Exception(f"{TTestCompared.LineCount} comparison type is not allowed for"
                                f" output type {output_type}")
        else:
            raise Exception(f"{expected_output.compared} comparison type is not supported")

        return value

    def read_md5(self, file_path: str):
        hash_md5 = hashlib.md5()
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
        return hash_md5.hexdigest()
