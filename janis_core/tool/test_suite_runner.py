import difflib
import hashlib
import os
import urllib
from typing import Dict, List, Any, Optional, Tuple, Set

from janis_core.tool.tool import Tool
from janis_core.tool.test_classes import (
    TTestExpectedOutput,
    TTestPreprocessor,
    TTestCase,
)
from janis_core.types import File, String, Array
from janis_core.tool import test_helpers
from janis_core.utils.secondary import apply_secondary_file_format_to_filename

from janis_core import Logger


class ToolTestSuiteRunner:
    """
    A class to read tool test cases, run the test cases, and assert the expected output
    """

    def __init__(self, tool: Tool, config: str = None):
        self.tool = tool
        self.output_dir = os.path.join(os.getcwd(), "tests_output", self.tool.id())
        self.cached_input_files_dir = os.path.join(
            os.getcwd(), "tests_output", "cached_test_files"
        )

        self._raw_config = config
        self._config = None

    @property
    def config(self):
        if self._config is None:

            from janis_assistant.management.configuration import JanisConfiguration

            if not self._raw_config:
                pass
            elif isinstance(self._raw_config, JanisConfiguration):
                self._config = self._raw_config
            elif isinstance(self._raw_config, str):
                self._config = JanisConfiguration.initial_configuration(
                    path=self._raw_config
                )
            else:
                raise ValueError(
                    f"Unrecognised type for janis configuration {self._raw_config} (type: {type(self._raw_config)})"
                )
        return self._config

    def run(self, input: Dict[str, str], engine: str) -> Dict[str, Any]:
        """
        Run a tool or workflow given a list of input using the specified engine

        :param input: workflow or tool input
        :type input: Dict[str, str]
        :param engine: name of engine to run this workflow/tool on
        :type engine: str

        :return:
        :rtype: Dict[str, Any]
        """
        # Note, we do not want janis_core to depend on janis_assistant except for when we run these tests
        # So, we only import this package here in this function
        # Make sure the correct version of janis_assistant is installed first
        test_helpers.verify_janis_assistant_installed()

        from janis_assistant.main import run_with_outputs

        output = run_with_outputs(
            tool=self.tool,
            inputs=input,
            output_dir=self.output_dir,
            engine=engine,
            config=self.config,
        )

        return output

    def run_one_test_case(
        self, t: TTestCase, engine: str, output: Optional[Dict] = None
    ) -> Tuple[Set, Set, Dict]:
        """
        Run one test case and assert multiple expected output

        :param t: test case
        :type t: TTestCase
        :param engine: name of engine to run this workflow/tool on
        :type engine: str

        :return: A tuple of failed error messages and successful expected output
        :rtype: Tuple[Set, Set]
        """
        # If output is provided, we don't want to re-run the workflow
        if output is None:
            output = self.run(input=t.input, engine=engine)

        failed = set()
        succeeded = set()

        for test_logic in t.output:
            workflow_output = output[test_logic.tag]
            self._download_remote_files(test_logic)
            actual_output = self.get_value_to_compare(test_logic, workflow_output)
            expected_value = self.get_expected_value(test_logic)
            test_result = test_logic.operator(actual_output, expected_value)

            if test_result is False:
                failed.add(
                    f"{str(test_logic)} {type(expected_value)}"
                    f" | actual output: {actual_output} {type(actual_output)}"
                )
            else:
                succeeded.add(str(test_logic))

        return failed, succeeded, output

    def get_expected_value(self, test_logic: TTestExpectedOutput) -> Any:
        """
        Get the expected value of a given expected output.
        It could be specified in a variable or in a file.

        :param test_logic: an object that holds information about an expected output
        :type test_logic: TTestExpectedOutput

        :return: expected value
        :rtype: Any
        """
        if test_logic.expected_value is not None:
            return test_logic.expected_value
        elif test_logic.expected_file is not None:
            with open(test_logic.expected_file) as expected_file:
                return expected_file.read()
        else:
            raise Exception(
                "one of `expected_value` or `expected_file` must not be empty"
            )

    def get_value_to_compare(
        self, test_logic: TTestExpectedOutput, output_value: Any
    ) -> Any:
        """
        convert workflow output value to the format we want to compare

        :param test_logic: an object that holds information about an expected output
        :type test_logic: TTestExpectedOutput
        :param output_value: output value from the executed workflow/tool
        :type output_value: Any

        :return:
        :rtype:
        """
        output_tag = test_logic.tag
        output_type = self.tool.outputs_map().get(output_tag).outtype
        output_value = self._extract_workflow_output(
            test_logic=test_logic, output_value=output_value, output_type=output_type
        )

        # Convert the output value to a format that we want to apply our test (e.g. md5, file content, etc)
        # Use a user-provided preprocessor or our predefined list of preprocessors
        if callable(test_logic.preprocessor):
            return test_logic.preprocessor(output_value)
        else:
            return self._apply_preprocessor(
                test_logic=test_logic,
                output_value=output_value,
                output_type=output_type,
            )

    def _apply_preprocessor(
        self, test_logic: TTestExpectedOutput, output_value: Any, output_type: Any
    ) -> Any:
        """
        Transform the output value based on the type of value we want to compare (e.g. md5, file size, etc)

        :param test_logic: an object that holds information about an expected output
        :type test_logic: TTestExpectedOutput
        :param output_value: output value from the executed workflow/tool
        :type output_value: Any
        :param output_type: Janis output data type
        :type output_type: Any

        :return: reformatted output value based on the specified preprocessor
        :rtype: Any
        """
        value = None

        if test_logic.preprocessor == TTestPreprocessor.Value:
            value = output_value
        elif test_logic.preprocessor == TTestPreprocessor.FileContent:
            with open(output_value) as f:
                value = f.read()
        elif test_logic.preprocessor == TTestPreprocessor.FileExists:
            value = os.path.isfile(output_value)
        elif test_logic.preprocessor == TTestPreprocessor.FileDiff:
            value = self.file_diff(
                expected_file_path=test_logic.file_diff_source,
                output_file_path=output_value,
            )
        elif test_logic.preprocessor == TTestPreprocessor.FileMd5:
            value = self.read_md5(output_value)
        elif test_logic.preprocessor == TTestPreprocessor.FileSize:
            value = os.path.getsize(output_value)
        elif test_logic.preprocessor == TTestPreprocessor.LineCount:
            value = self.line_count(output_type=output_type, output_value=output_value)
        elif test_logic.preprocessor == TTestPreprocessor.ListSize:
            if not output_value:
                return 0
            value = len(output_value.split("|"))
        else:
            raise Exception(
                f"{test_logic.preprocessor} comparison type is not supported"
            )

        return value

    def _extract_workflow_output(
        self, test_logic: TTestExpectedOutput, output_value: Any, output_type: Any
    ) -> Any:
        """
        When we want to apply test to other output related to this output.

        e.g. converting from a list of values to a single value,
        Append suffix to filenames if we want to check secondary files, etc

        :param test_logic: an object that holds information about an expected output
        :type test_logic: TTestExpectedOutput
        :param output_value: output value from the executed workflow/tool
        :type output_value: Any
        :param output_type: Janis output data type
        :type output_type: Any

        :return: value extracted from workflow/tool output based on the defined test logic
        :rtype: Any
        """
        # Convert array to element of array
        if isinstance(output_type, Array):
            if test_logic.array_index is not None:
                output_list = output_value.split("|")
                output_value = output_list[test_logic.array_index]

        # Add extension to a filename (when testing secondary files)
        if isinstance(output_type, File):
            if test_logic.suffix is not None:
                output_value = apply_secondary_file_format_to_filename(
                    output_value, test_logic.suffix
                )

        return output_value

    def _download_remote_files(self, test_logic: TTestExpectedOutput):
        """
        Download remote test files (only expected output files) to a cache directory

        :param test_logic: an object that holds information about an expected output
        :type test_logic: TTestExpectedOutput
        :return: None
        :rtype: None
        """

        file_attributes = ["expected_file", "file_diff_source"]
        for att in file_attributes:
            if not hasattr(test_logic, att):
                raise Exception(f"{test_logic.__class__} has no attribute {att}")

            source = getattr(test_logic, att)

            if source:

                test_helpers.verify_janis_assistant_installed()
                from janis_assistant.management.filescheme import (
                    FileScheme,
                    LocalFileScheme,
                )

                # f = FileScheme(source)
                if not FileScheme.is_local_path(source):
                    fs = FileScheme.get_filescheme_for_url(source)
                    last_modified = fs.last_modified(source)

                    local_file_path = os.path.join(
                        self.cached_input_files_dir,
                        f"{test_helpers.hash_filename(source, last_modified)}_{os.path.basename(source)}",
                    )

                    # Only download if the file does not already exist
                    if not os.path.exists(local_file_path):
                        Logger.info(f"Downloading remote file to {local_file_path}")

                        os.makedirs(self.cached_input_files_dir, exist_ok=True)
                        fs.cp_from(source, local_file_path)
                    else:
                        Logger.info(
                            f"Skip downloading remote file. File {source} already exists in {local_file_path}"
                        )

                setattr(test_logic, att, local_file_path)

    def read_md5(self, file_path: str) -> str:
        """
        Read md5 checksum of a file

        :param file_path: full path to a file
        :type file_path: str

        :return: md5 checksum of a file
        :rtype: str
        """
        hash_md5 = hashlib.md5()
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
        return hash_md5.hexdigest()

    def file_diff(
        self,
        expected_file_path: str,
        output_file_path: Optional[str] = None,
        output_content: Optional[str] = None,
    ) -> List[str]:
        """
        Create a list if file diff result between two files

        :param expected_file_path: path to the source file to perform diff with
        :type expected_file_path: str
        :param output_file_path: path to the workflow/tool output file to perform diff with
        :type output_file_path: Optional[str]
        :param output_content: a string from workflow/tool output to perform the diff with
        :type output_content: Optional[str]

        :return: a list of diff output
        :rtype: List[str]
        """
        with open(expected_file_path) as expected_file:
            expected_content = list(expected_file)

        if output_file_path is not None:
            with open(output_file_path) as output_file:
                output_content = list(output_file)

        if output_content is None:
            raise Exception(
                "One of `output_file_path` or `output_content` must be provided"
            )

        diff = difflib.unified_diff(
            expected_content,
            output_content,
            fromfile="expected",
            tofile="actual",
            lineterm="",
        )
        return list(diff)

    def line_count(self, output_value: Any, output_type: Any) -> int:
        """
        Count number of lines in a file or in a string

        :param output_value: workflow/tool output value
        :type output_value: Any
        :param output_type: Janis data type
        :type output_type: Any

        :return: number of lines
        :rtype: int
        """
        if isinstance(output_type, File):
            # text file only here
            with open(output_value) as fp:
                value = sum(1 for _ in fp)
        elif isinstance(output_type, String):
            value = len(output_value.splitlines())
        else:
            raise Exception(
                f"{TTestPreprocessor.LineCount} comparison type is not allowed for"
                f" output type {output_type}"
            )

        return value
