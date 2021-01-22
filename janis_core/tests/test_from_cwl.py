import unittest

from janis_core import CWlParser


class TestFromCwlExpressions(unittest.TestCase):
    def test_number(self):
        result = CWlParser(cwl_version="v1.0").parse_basic_expression("${return 16 }")
        self.assertEqual(16, result)

    def test_number_with_semicolon(self):
        result = CWlParser(cwl_version="v1.0").parse_basic_expression("${return 16;}")
        self.assertEqual(16, result)

    def test_number_with_spaces(self):
        result = CWlParser(cwl_version="v1.0").parse_basic_expression(
            "${ return 80000 }"
        )
        self.assertEqual(80000, result)

    def test_string(self):
        result = CWlParser(cwl_version="v1.0").parse_basic_expression(
            '${ return "over 80000" }'
        )
        self.assertEqual("over 80000", result)
