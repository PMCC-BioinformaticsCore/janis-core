import difflib
import hashlib
import os
from typing import Dict, List, Any

from janis_core.tool.tool import Tool
from janis_core.tool.test_classes import TTestExpectedOutput, TTestCompared, TTestCase
from janis_core.types import File, String, Array
from janis_core.tool import test_helpers


class ToolTestSuiteRunner:
    """
    A class to read tool test cases, run the test cases, and assert the expected output
    """

    def __init__(self, tool: Tool):
        self.tool = tool

    def run(self, input: Dict[str, str], engine):
        # Note, we do not want janis_core to depend on janis_assistant except for when we run these tests
        # So, we only import this package here in this function
        # Make sure the correct version of janis_assistant is installed first

        test_helpers.verify_janis_assistant_installed()

        from janis_assistant.main import run_with_outputs
        from janis_assistant.management.configuration import JanisConfiguration

        output_dir = os.path.join(os.getcwd(), "tests_output", self.tool.id())
        config = JanisConfiguration(engine=engine)
        output = run_with_outputs(
            tool=self.tool, inputs=input, output_dir=output_dir, config=config
        )

        return output

    def run_one_test_case(self, t: TTestCase, engine):
        """
        Run one test case and assert multiple expected output
        """
        output = self.run(input=t.input, engine=engine)

        failed = set()
        succeeded = set()

        for test_logic in t.output:
            workflow_output = output[test_logic.tag]
            actual_output = self.get_value_to_compare(test_logic, workflow_output)
            test_result = test_logic.operator(actual_output, test_logic.expected_value)

            if test_result is False:
                failed.add(
                    f"{str(test_logic)} {type(test_logic.expected_value)}"
                    f" | actual output: {actual_output} {type(actual_output)}"
                )
            else:
                succeeded.add(str(test_logic))

        return failed, succeeded

    def get_value_to_compare(
        self, expected_output: TTestExpectedOutput, output_value: Any
    ) -> Any:
        """
        convert workflow output value to the format we want to compare
        """
        output_tag = expected_output.tag
        output_type = self.tool.outputs_map().get(output_tag).outtype

        # Convert array to element of array
        if isinstance(output_type, Array):
            if expected_output.array_index is not None:
                output_list = output_value.split("|")
                output_value = output_list[expected_output.array_index]

        # Add extension to a filename (when testing secondary files)
        if isinstance(output_type, File):
            if expected_output.suffix is not None:
                output_value = f"{output_value}{expected_output.suffix}"

        value = None

        if expected_output.compared == TTestCompared.Value:
            value = output_value
        elif expected_output.compared == TTestCompared.FileContent:
            with open(output_value) as f:
                value = f.read()
        elif expected_output.compared == TTestCompared.FileExists:
            value = os.path.isfile(output_value)
        elif expected_output.compared == TTestCompared.FileDiff:
            value = self.file_diff(
                expected_file_path=expected_output.expected_source,
                output_file_path=output_value,
            )
        elif expected_output.compared == TTestCompared.FileMd5:
            value = self.read_md5(output_value)
        elif expected_output.compared == TTestCompared.FileSize:
            value = os.path.getsize(output_value)
        elif expected_output.compared == TTestCompared.LineCount:
            value = self.line_count(output_type=output_type, output_value=output_value)
        elif expected_output.compared == TTestCompared.ListSize:
            pass
        else:
            raise Exception(
                f"{expected_output.compared} comparison type is not supported"
            )

        return value

    def read_md5(self, file_path: str) -> str:
        hash_md5 = hashlib.md5()
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
        return hash_md5.hexdigest()

    def file_diff(self, expected_file_path: str, output_file_path: str) -> List[str]:
        with open(expected_file_path) as expected_file:
            expected_content = list(expected_file)
        with open(output_file_path) as output_file:
            output_content = list(output_file)

        diff = difflib.unified_diff(
            expected_content,
            output_content,
            fromfile="expected",
            tofile="actual",
            lineterm="",
        )
        return list(diff)

    def line_count(self, output_value: Any, output_type: Any) -> int:
        try:
            if isinstance(output_type, File):
                # text file only here
                value = sum(1 for line in open(output_value, "r"))
            elif isinstance(output_type, String):
                value = len(output_value.splitlines())
        except Exception as e:
            raise Exception(
                f"{TTestCompared.LineCount} comparison type is not allowed for"
                f" output type {output_type}"
            )

        return value
