from enum import Enum
from typing import Optional, Dict, Any, Callable


class TTestCompared(Enum):
    Value = "value"
    FileDiff = "file-diff"
    FileContent = "file-content"
    FileExists = "file-exists"
    FileSize = "file-size"
    FileMd5 = "file-md5"
    LineCount = "line-count"
    ListSize = "list-size"


class TTestExpectedOutput(object):
    def __init__(
        self,
        tag: str,
        compared: Enum,
        operator: Callable,
        expected_value: Optional[Any] = None,
        expected_file: Optional[str] = None,
        file_diff_source: Optional[str] = None,
        array_index: Optional[int] = None,
        suffix: Optional[str] = None,
    ):
        self.tag = tag
        self.compared = compared
        self.operator = operator
        self.expected_value = expected_value
        self.expected_file = expected_file

        # If the 'compared' type requires us to transform from a different object to expected value
        # example:
        # compared: FileDiff
        # file_diff_source: file path to the expected file
        # expected_value: values readable by 'operator'. It could be the number of
        # lines of diff result of expected file and actual file.
        self.file_diff_source = file_diff_source

        # if an output is an array, we can look at just 1 item of the array, so here we specify the index
        self.array_index = array_index

        # if the compared object is a file, we can add suffix to test secondary files of this file
        self.suffix = suffix

    def __repr__(self):
        repr_expected_value = str(self.expected_value)
        if self.expected_value is None:
            repr_expected_value = f"content of {self.expected_file}"

        return f"{self.tag}: {self.compared.value} {str(self.operator)} {repr_expected_value}"


class TTestCase(object):
    def __init__(
        self, name: str, input: Dict[str, Any], output: Dict[str, TTestExpectedOutput]
    ):
        self.name = name
        self.input = input
        self.output = output
