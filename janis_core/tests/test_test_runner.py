import os
import operator
from typing import Optional, List, Dict, Union
from unittest import TestCase, mock
from janis_core.tool.test_suite_runner import ToolTestSuiteRunner
from janis_core.tool.test_classes import (
    TTestCase,
    TTestExpectedOutput,
    TTestPreprocessor,
)


from janis_core import ToolOutput, ToolInput, CommandTool, Stdout, File

valid_url = "https://abc.com/some_dir/expected_output_file.txt"
valid_url_2 = "https://abc.com/some_dir/diff_file.txt"
valid_url_last_modified = "Thu, 10 Dec 2020 02:07:45 GMT"
valid_url_content = "test data content"
valid_url_content_2 = "test diff content"
expected_tool_output = "some expected output text"
expected_file = "/my/local/file"


class MockResponse:
    def __init__(
        self,
        json_data: Dict,
        status_code: int,
        headers: Optional[Dict] = None,
        content: Optional[str] = None,
    ):
        self.json_data = json_data
        self.status_code = status_code
        self.headers = headers
        self.content = content.encode()

    def json(self):
        return self.json_data


# This method will be used by the mock to replace requests.get
def mocked_requests_head_and_get(*args, **kwargs):

    if args[0] == valid_url:
        return MockResponse(
            {}, 200, {"Last-Modified": valid_url_last_modified}, valid_url_content
        )
    elif args[0] == valid_url_2:
        return MockResponse(
            {}, 200, {"Last-Modified": valid_url_last_modified}, valid_url_content_2
        )

    return MockResponse(None, 404)


def mocked_run_with_outputs(*args, **kwargs):
    return {"tool_output": expected_tool_output, "file": expected_file}


class TestTool(CommandTool):
    def tool(self) -> str:
        return "EchoTestTool"

    def base_command(self) -> Optional[Union[str, List[str]]]:
        return "echo"

    def inputs(self) -> List[ToolInput]:
        return [ToolInput("inp", str, position=0)]

    def outputs(self):
        return [ToolOutput("tool_output", Stdout), ToolOutput("file", Stdout)]

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
            tag="tool_output",
            preprocessor=TTestPreprocessor.Value,
            operator=operator.eq,
            expected_value=5,
        )
        assert runner.get_expected_value(t1) == 5

        t2 = TTestExpectedOutput(
            tag="tool_output",
            preprocessor=TTestPreprocessor.Value,
            operator=operator.eq,
            expected_value=["abc", "def"],
        )
        assert runner.get_expected_value(t2) == ["abc", "def"]

    def test_get_value_to_compare(self):
        runner = ToolTestSuiteRunner(self.tool)

        t1 = TTestExpectedOutput(
            tag="tool_output",
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

    @mock.patch("requests.get", side_effect=mocked_requests_head_and_get)
    @mock.patch("requests.head", side_effect=mocked_requests_head_and_get)
    def test_download_remote_files(self, mock_get, mock_head):
        runner = ToolTestSuiteRunner(self.tool)

        t1 = TTestExpectedOutput(
            tag="tool_output",
            preprocessor=TTestPreprocessor.Value,
            operator=operator.eq,
            expected_file=valid_url,
            file_diff_source=valid_url_2,
        )

        f1_local = f"f282f7ae4b77890d7c4f27740d1507e1_{os.path.basename(valid_url)}"
        f1_path = os.path.join(
            os.getcwd(),
            "tests_output",
            "cached_test_files",
            f1_local,
        )

        f2_local = f"c8ea947c092f9b75da21d24420c8a3a1_{os.path.basename(valid_url_2)}"
        f2_path = os.path.join(
            os.getcwd(),
            "tests_output",
            "cached_test_files",
            f2_local,
        )

        assert t1.expected_file == valid_url
        assert t1.file_diff_source == valid_url_2
        runner._download_remote_files(t1)

        assert t1.expected_file == f1_path
        assert os.path.exists(f1_path)

        with open(f1_path, "rb") as f:
            content = f.read()

        assert content == valid_url_content.encode()

        assert t1.file_diff_source == f2_path
        assert os.path.exists(f2_path)

        with open(f2_path, "rb") as f:
            content_2 = f.read()

        assert content_2 == valid_url_content_2.encode()

    @mock.patch(
        "janis_assistant.main.run_with_outputs",
        side_effect=mocked_run_with_outputs,
    )
    def test_run_one_test_case_succeed(self, mock_output):
        runner = ToolTestSuiteRunner(self.tool)

        tc = TTestCase(
            name="test case 1",
            input={"inp": "my input"},
            output=[
                TTestExpectedOutput(
                    tag="tool_output",
                    preprocessor=TTestPreprocessor.Value,
                    operator=operator.eq,
                    expected_value=expected_tool_output,
                )
            ],
        )

        failed, succeeded, output = runner.run_one_test_case(tc, engine="cromwell")

        assert output == mock_output()
        assert failed == set()
        assert succeeded == {"tool_output: value eq some expected output text"}

    @mock.patch(
        "janis_assistant.main.run_with_outputs",
        side_effect=mocked_run_with_outputs,
    )
    def test_run_one_test_case_fail(self, mock_output):
        runner = ToolTestSuiteRunner(self.tool)

        tc = TTestCase(
            name="test case 1",
            input={"inp": "my input"},
            output=[
                TTestExpectedOutput(
                    tag="tool_output",
                    preprocessor=TTestPreprocessor.Value,
                    operator=operator.eq,
                    expected_value=expected_tool_output,
                ),
                TTestExpectedOutput(
                    tag="file",
                    preprocessor=TTestPreprocessor.Value,
                    operator=operator.eq,
                    expected_value="/xxx/yyy",
                ),
            ],
        )

        failed, succeeded, output = runner.run_one_test_case(tc, engine="cromwell")

        assert output == mock_output()
        assert failed == {
            "file: value eq /xxx/yyy <class 'str'> | actual output: /my/local/file <class 'str'>"
        }
        assert succeeded == {"tool_output: value eq some expected output text"}

    @mock.patch(
        "janis_assistant.main.run_with_outputs",
        side_effect=mocked_run_with_outputs,
    )
    def test_run_one_test_case_dry_run(self, mock_output):
        runner = ToolTestSuiteRunner(self.tool)

        tc = TTestCase(
            name="test case 1",
            input={"inp": "my input"},
            output=[
                TTestExpectedOutput(
                    tag="tool_output",
                    preprocessor=TTestPreprocessor.Value,
                    operator=operator.eq,
                    expected_value=expected_tool_output,
                )
            ],
        )

        expected_output = {"tool_output": "from dry run"}
        failed, succeeded, output = runner.run_one_test_case(
            tc, engine="cromwell", output=expected_output
        )
        print(failed)

        assert output == expected_output
        assert failed == {
            "tool_output: value eq some expected output text <class 'str'> | actual output: from dry run <class 'str'>"
        }
        assert succeeded == set()
