import unittest

from janis_core import Array, String, Int


class TestValidateInputvalue(unittest.TestCase):
    def test_validate_string_optional_allowoptional_value(self):
        self.assertTrue(String(optional=True).validate_value("aa", True))

    def test_validate_string_optional_allowoptional_novalue(self):
        self.assertTrue(String(optional=True).validate_value(None, True))

    def test_validate_string_optional_disallowoptional_value(self):
        self.assertTrue(String(optional=True).validate_value("aa", False))

    def test_validate_string_optional_disallowoptional_novalue(self):
        self.assertTrue(String(optional=True).validate_value(None, False))

    def test_validate_string_nooptional_allowoptional_value(self):
        self.assertTrue(String().validate_value("aa", True))

    def test_validate_string_nooptional_allowoptional_novalue(self):
        self.assertTrue(String().validate_value(None, True))

    def test_validate_string_nooptional_disallowoptional_value(self):
        self.assertTrue(String().validate_value("aa", False))

    def test_validate_string_nooptional_disallowoptional_novalue(self):
        self.assertFalse(String().validate_value(None, False))

    def test_array_valid(self):
        self.assertTrue(Array(String()).validate_value(["aa", "bb"], True))

    def test_array_valid_optional_internal(self):
        self.assertTrue(
            Array(String(optional=True)).validate_value(["aa", None], False)
        )

    def test_array_invalid_int_string(self):
        self.assertTrue(Array(String()).validate_value(["aa", 2], True))

    def test_int_invalid_float(self):
        self.assertFalse(Array(Int()).validate_value(2.0, True))
