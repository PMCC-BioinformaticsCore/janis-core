import operator
import os
from typing import Optional, List, Dict, Union
from unittest import TestCase, mock
from janis_core.tool.test_suite_runner import ToolTestSuiteRunner
from janis_core.tool.test_classes import (
    TTestCase,
    TTestExpectedOutput,
    TTestPreprocessor,
)
from janis_core.types import File, String, Array

from janis_core import ToolOutput, ToolInput, CommandTool, Stdout
from janis_core.tool.test_classes import TTestExpectedOutput, TTestPreprocessor
from janis_core.tool.test_suite_runner import ToolTestSuiteRunner
from nose.tools import nottest


valid_url = "https://abc.com/some_dir/expected_output_file.txt"
valid_url_2 = "https://abc.com/some_dir/diff_file.txt"
valid_url_last_modified = "Thu, 10 Dec 2020 02:07:45 GMT"
valid_url_content = "test data content"
valid_url_content_2 = "test diff content"
expected_tool_output = "some expected output text"
expected_file = "/my/local/file"


class MockResponse:
    def __init__(
        self, json_data: Dict, status_code: int, headers: Optional[Dict] = None
    ):
        self.json_data = json_data
        self.status_code = status_code
        self.headers = headers

    def getcode(self):
        return self.status_code

    def getheader(self, field: str):
        return self.headers.get(field)

    def json(self):
        return self.json_data


def mocked_urllib_urlopen(*args, **kwargs):
    if args[0].full_url == valid_url:
        return MockResponse({}, 200, {"Last-Modified": valid_url_last_modified})
    elif args[0].full_url == valid_url_2:
        return MockResponse({}, 200, {"Last-Modified": valid_url_last_modified})

    return MockResponse(None, 404)


def mocked_urllib_urlretrieve(*args, **kwargs):

    if kwargs["url"] == valid_url:
        local_path = kwargs["filename"]
        with open(local_path, "wb") as f:
            f.write(valid_url_content.encode())
        return (local_path, {"Last-Modified": valid_url_last_modified})

    elif kwargs["url"] == valid_url_2:
        local_path = kwargs["filename"]
        with open(local_path, "wb") as f:
            f.write(valid_url_content_2.encode())
        return (local_path, {"Last-Modified": valid_url_last_modified})

    return None, None


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
    """
    This class is just testing the 'expected_output' evaluator, and not actually running tools.
    If you're interested in running tools (and require janis_assistant, consider the following)
        from unittest import TestCase, skipUnless
        try:
            import janis_assistant

            has_janis_assistant = True
        except:
            has_janis_assistant = False

        # ...
        class MyTest(TestCase):
            @skipUnless(has_janis_assistant, "Must have janis assistant to run this test")
            def test_thing_one(self):
                self.assertTrue
    """

    def setUp(self):
        self.tool = TestTool()
        self.test_data_dir = os.path.join(
            os.path.dirname(os.path.realpath(__file__)), "data"
        )

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

    @nottest
    @mock.patch("urllib.request.urlopen", side_effect=mocked_urllib_urlopen)
    @mock.patch("urllib.request.urlretrieve", side_effect=mocked_urllib_urlretrieve)
    def test_download_remote_files(self, mock_urlopen, mock_urlretrieve):
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

    def test_run_one_test_case_succeed(
        self,
    ):
        runner = ToolTestSuiteRunner(self.tool)

        tc = TTestCase(
            name="test1",
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

        expected_output = {"tool_output": expected_tool_output, "file": expected_file}
        failed, succeeded, output = runner.run_one_test_case(
            tc, engine="cromwell", output=expected_output
        )

        assert output == expected_output
        assert failed == []
        assert succeeded == ["tool_output: value eq some expected output text"]

    def test_run_one_test_case_fail(self):
        runner = ToolTestSuiteRunner(self.tool)

        tc = TTestCase(
            name="test1",
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

        expected_output = {"tool_output": expected_tool_output, "file": expected_file}
        failed, succeeded, output = runner.run_one_test_case(
            tc, engine="cromwell", output=expected_output
        )

        assert output == expected_output
        assert failed == [
            "file: value eq /xxx/yyy <class 'str'> | actual output: /my/local/file <class 'str'>"
        ]
        assert succeeded == ["tool_output: value eq some expected output text"]

    def test_apply_preprocessor(self):
        runner = ToolTestSuiteRunner(self.tool)

        t = TTestExpectedOutput(
            tag="tool_output",
            preprocessor=TTestPreprocessor.Value,
            operator=operator.eq,
            expected_value="xxx",
        )

        assert runner._apply_preprocessor(t, "aaa", str) == "aaa"

        file_path = os.path.join(self.test_data_dir, "test.txt")
        t = TTestExpectedOutput(
            tag="tool_output",
            preprocessor=TTestPreprocessor.FileContent,
            operator=operator.eq,
            expected_value="xxx",
        )

        assert runner._apply_preprocessor(t, file_path, File()) == "abc\ndef\nghi\n"

        file_path = os.path.join(self.test_data_dir, "test.txt")
        t = TTestExpectedOutput(
            tag="tool_output",
            preprocessor=TTestPreprocessor.FileExists,
            operator=operator.eq,
            expected_value="xxx",
        )

        assert runner._apply_preprocessor(t, file_path, File()) is True

        expected_file_path = os.path.join(self.test_data_dir, "test.txt")
        output_file_path = os.path.join(self.test_data_dir, "test_diff.txt")
        t = TTestExpectedOutput(
            tag="tool_output",
            preprocessor=TTestPreprocessor.FileDiff,
            operator=operator.eq,
            file_diff_source=expected_file_path,
            expected_value="",
        )

        assert runner._apply_preprocessor(t, output_file_path, File()) == [
            "--- expected",
            "+++ actual",
            "@@ -1,3 +1,2 @@",
            " abc\n",
            "-def\n",
            " ghi\n",
        ]

        file_path = os.path.join(self.test_data_dir, "test.txt")
        t = TTestExpectedOutput(
            tag="tool_output",
            preprocessor=TTestPreprocessor.FileMd5,
            operator=operator.eq,
            expected_value="xxx",
        )

        assert (
            runner._apply_preprocessor(t, file_path, File())
            == "7fe4b98fb68caf83fa31905435f08da0"
        )

        file_path = os.path.join(self.test_data_dir, "test.txt")
        t = TTestExpectedOutput(
            tag="tool_output",
            preprocessor=TTestPreprocessor.FileSize,
            operator=operator.eq,
            expected_value="xxx",
        )

        assert runner._apply_preprocessor(t, file_path, File()) == 12

        file_path = os.path.join(self.test_data_dir, "test.txt")
        t = TTestExpectedOutput(
            tag="tool_output",
            preprocessor=TTestPreprocessor.LineCount,
            operator=operator.eq,
            expected_value="xxx",
        )

        assert runner._apply_preprocessor(t, file_path, File()) == 3

        t = TTestExpectedOutput(
            tag="tool_output",
            preprocessor=TTestPreprocessor.LineCount,
            operator=operator.eq,
            expected_value="xxx",
        )

        assert runner._apply_preprocessor(t, "line1\nline2", String()) == 2
        assert runner._apply_preprocessor(t, "line1\nline2\n", String()) == 2
        assert runner._apply_preprocessor(t, "line1", String()) == 1
        assert runner._apply_preprocessor(t, "", String()) == 0

        t = TTestExpectedOutput(
            tag="tool_output",
            preprocessor=TTestPreprocessor.ListSize,
            operator=operator.eq,
            expected_value="xxx",
        )

        assert runner._apply_preprocessor(t, "a|b|c|d|e", Array(String)) == 5
        assert runner._apply_preprocessor(t, "", Array(String)) == 0

    def test_extract_workflow_output(self):
        runner = ToolTestSuiteRunner(self.tool)

        t = TTestExpectedOutput(
            tag="tool_output",
            preprocessor=TTestPreprocessor.Value,
            operator=operator.eq,
            expected_value="xxx",
            array_index=2,
        )

        assert runner._extract_workflow_output(t, "a|b|c|d|e", Array(String)) == "c"

        t = TTestExpectedOutput(
            tag="tool_output",
            preprocessor=TTestPreprocessor.Value,
            operator=operator.eq,
            expected_value="xxx",
            suffix_secondary_file=".bai",
        )

        assert (
            runner._extract_workflow_output(t, "/dir/sample.bam", File())
            == "/dir/sample.bam.bai"
        )

        t = TTestExpectedOutput(
            tag="tool_output",
            preprocessor=TTestPreprocessor.Value,
            operator=operator.eq,
            expected_value="xxx",
            suffix_secondary_file="^.bai",
        )

        assert (
            runner._extract_workflow_output(t, "/dir/sample.bam", File())
            == "/dir/sample.bai"
        )
