import unittest

from janis_core import (
    CWlParser,
    InputSelector,
    BasenameOperator,
    FileSizeOperator,
    ReadContents,
)


class TestFromCwlExpressions(unittest.TestCase):
    parser = CWlParser(cwl_version="v1.0")

    def test_number(self):
        expr = "${return 16 }"
        result = self.parser.parse_basic_expression(expr)
        self.assertEqual(16, result)

    def test_number_with_semicolon(self):
        expr = "${return 16;}"
        result = self.parser.parse_basic_expression(expr)
        self.assertEqual(16, result)

    def test_number_with_spaces(self):
        expr = "${ return 80000 }"
        result = self.parser.parse_basic_expression(expr)
        self.assertEqual(80000, result)

    def test_string(self):
        expr = '${ return "over 80000" }'
        result = self.parser.parse_basic_expression(expr)
        self.assertEqual("over 80000", result)

    def test_input_selector(self):
        expr = "$(inputs.my_input)"
        result = self.parser.parse_basic_expression(expr)
        self.assertIsInstance(result, InputSelector)
        self.assertEqual("my_input", result.input_to_select)

    def test_input_selector_with_basename(self):
        expr = "$(inputs.my_input.basename)"
        result = self.parser.parse_basic_expression(expr)
        self.assertIsInstance(result, BasenameOperator)
        self.assertIsInstance(result.args[0], InputSelector)
        self.assertEqual("my_input", result.args[0].input_to_select)

    def test_input_selector_with_filesize(self):
        expr = "$(inputs.my_input.size)"
        result = self.parser.parse_basic_expression(expr)
        self.assertIsInstance(result, FileSizeOperator)
        self.assertIsInstance(result.args[0], InputSelector)
        self.assertEqual("my_input", result.args[0].input_to_select)

    def test_input_selector_with_contents(self):
        expr = "$(inputs.my_input.contents)"
        result = self.parser.parse_basic_expression(expr)
        self.assertIsInstance(result, ReadContents)
        self.assertIsInstance(result.args[0], InputSelector)
        self.assertEqual("my_input", result.args[0].input_to_select)
