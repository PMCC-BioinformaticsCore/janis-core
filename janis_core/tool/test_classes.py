import re
from enum import Enum
from typing import Optional, Dict, Any, Callable, Union, List


class TTestPreprocessor(Enum):
    """
    An Enum class for a list of possible pre-processing we can do to output values
    """

    Value = "value"
    FileDiff = "file-diff"
    FileContent = "file-content"
    FileExists = "file-exists"
    FileSize = "file-size"
    FileMd5 = "file-md5"
    LineCount = "line-count"
    ListSize = "list-size"


class TTestExpectedOutput(object):
    """
    Describe the logic on how to test the expected output of a test case.
    A test case can have multiple instances of this class to test different output or
    different logic of the same output
    """

    def __init__(
        self,
        tag: str,
        preprocessor: Union[TTestPreprocessor, Callable[[Any], Any]],
        operator: Callable[[Any, Any], bool],
        expected_value: Optional[Any] = None,
        expected_file: Optional[str] = None,
        file_diff_source: Optional[str] = None,
        array_index: Optional[int] = None,
        suffix_secondary_file: Optional[str] = None,
    ):
        """

        :param tag: output field name
        :type tag: str
        :param preprocessor: one of the TTestValuePreProcessor or a function to do the pre-processing of the output value
        :type preprocessor: Union[TTestPreprocessor, Callable]
        :param operator: A callable function to compare output value and expected value
        :type operator: Callable[[Any, Any], bool]
        :param expected_value: the value expected output
        :type expected_value:Optional[Any]
        :param expected_file: the file path to the value of expected output
        :type expected_file: Optional[str]
        :param file_diff_source: the file path to the source file for file diff comparison
        :type file_diff_source: Optional[str]
        :param array_index: an integer to represent the index of the element we want to test (for array output only)
        :type array_index: Optional[int]
        :param suffix_secondary_file: additional file
        :type suffix_secondary_file: Optional[str]
        """
        self.tag = tag
        self.preprocessor = preprocessor
        self.operator = operator
        self.expected_value = expected_value

        # we can have a 'string' expected output stored in a file
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
        self.suffix = suffix_secondary_file

        self._validate_input()

    def __repr__(self):
        repr_expected_value = str(self.expected_value)
        if self.expected_value is None:
            repr_expected_value = f"content of {self.expected_file}"

        repr_preprocessor = str(self.preprocessor)
        if isinstance(self.preprocessor, TTestPreprocessor):
            repr_preprocessor = self.preprocessor.value

        return f"{self.tag}: {repr_preprocessor} {str(self.operator.__name__)} {repr_expected_value}"

    def _validate_input(self):
        if self.expected_value is None and self.expected_file is None:
            raise Exception(
                "one of `expected_value` or `expected_file` must not be empty"
            )


class TTestCase(object):
    """
    A test case requires a workflow or tool to be run once (per engine).
    But, we can have multiple output to apply different test logic.
    """

    def __init__(
        self, name: str, input: Dict[str, Any], output: List[TTestExpectedOutput]
    ):
        """

        :param name: name of the test case
        :type name: str
        :param input: tool input values keyed by its input field name
        :type input: Dict[str, Any]
        :param output: List of definitions on how to assert the expected output
        :type output: List[TTestExpectedOutput]
        """
        # we don't want white space in the test case name
        if bool(re.search("\s", name)):
            raise Exception(f"No whitespace is allowed in the TestCase name {name}")

        self.name = name
        self.input = input
        self.output = output
