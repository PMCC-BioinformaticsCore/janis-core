import difflib
import hashlib
import os
from typing import Dict, List, Optional, Any
from pkg_resources import parse_version

from janis_core.tool.tool import Tool, TTestExpectedOutput, TTestCompared, TTestCase
from janis_core.types import DataType, File, String
from janis_core.tool import test_helpers


class ToolTestSuiteRunner():
    """
    A class to read tool test cases, run the test cases, and assert the expected output
    """
    def __init__(self, tool: Tool):
        self.tool = tool

    def run(self, input: Dict[str, str]):
        # Note, we do not want janis_core to depend on janis_assistant except for when we run these tests
        # So, we only import this package here in this function
        # Make sure the correct version of janis_assistant is installed first
        min_version_required = test_helpers.janis_assistant_version_required_min
        try:
            import janis_assistant
            if parse_version(janis_assistant.__version__) < parse_version(min_version_required):
                raise Exception(f"to run this test, janis_asisstant >= {min_version_required}"
                                f" must be installed. Installed version is {janis_assistant.__version__}")

            from janis_assistant.main import run_with_outputs
            from janis_assistant.management.configuration import JanisConfiguration
            from janis_assistant.engines.enginetypes import EngineType

            engines = [
                EngineType.cromwell.value,
                # EngineType.cwltool.value,
            ]
            current_dir = os.getcwd()

            output = {}
            for engine in engines:
                output_dir = os.path.join(current_dir, "tests_output", self.tool.id())
                config = JanisConfiguration(engine=engine)
                result = run_with_outputs(tool=self.tool, inputs=input, output_dir=output_dir, config=config)
                print(result)
                output[engine] = result

                #TODO: investigate why we need to change directory here?
                os.chdir(current_dir)

            return output
        except Exception as e:
            raise e

    def run_one_test_case(self, t: TTestCase):
        """
        Run one test case and assert multiple expected output
        """
        output_all_engines = self.run(t.input)

        test_results = {}
        for engine in output_all_engines:
            output = output_all_engines[engine]

            failed = set()
            succeeded = set()
            for expected_output in t.output:
                print(output)
                output_value = output[expected_output.tag]
                actual_output = self.get_value_to_compare(expected_output, output_value)

                test_result = expected_output.operator(actual_output, expected_output.expected_value)
                if test_result is False:
                    failed.add(f"{str(expected_output)} {type(expected_output.expected_value)} | actual output: {actual_output} {type(actual_output)}")
                else:
                    succeeded.add(str(expected_output))

            test_results[engine] = (failed, succeeded)

        return test_results

    def get_value_to_compare(self, expected_output: TTestExpectedOutput, output_value: Any) -> Any:
        output_tag = expected_output.tag
        value = None

        if expected_output.compared == TTestCompared.Value:
            value = output_value
        elif expected_output.compared == TTestCompared.FileContent:
            with open(output_value) as f:
                value = f.read()
        elif expected_output.compared == TTestCompared.FileDiff:
            value = self.file_diff(expected_file_path=expected_output.expected_source,
                                   output_file_path=output_value)
        elif expected_output.compared == TTestCompared.FileMd5:
            value = self.read_md5(output_value)
        elif expected_output.compared == TTestCompared.FileSize:
            value = os.path.getsize(output_value)
        elif expected_output.compared == TTestCompared.LineCount:
            value = self.line_count(output_tag=output_tag, output_value=output_value)
        else:
            raise Exception(f"{expected_output.compared} comparison type is not supported")

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

        diff = difflib.unified_diff(expected_content, output_content,
                                    fromfile='expected', tofile='actual', lineterm='')
        return list(diff)

    def line_count(self, output_tag: str, output_value: Any) -> int:
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

        return value
