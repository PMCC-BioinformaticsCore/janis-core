import unittest

from janis_core import Array, String, Stdout, File, Int, Float, Boolean
from janis_core.types import get_instantiated_type


class DataTypeWithSecondary(File):
    @staticmethod
    def name() -> str:
        return "test_secondary"

    @staticmethod
    def secondary_files():
        return [".txt"]


class TestTypes(unittest.TestCase):
    def test_array_of_strings(self):

        ar = Array(String())
        d = ar.cwl_type()
        self.assertEqual(d.get_dict(), {"type": "array", "items": "string"})

    def test_array_of_array_of_strings(self):
        ar = Array(Array(String()))
        d = ar.cwl_type()
        self.assertEqual(
            d.get_dict(),
            {"type": "array", "items": {"type": "array", "items": "string"}},
        )

    def test_stdout_normal(self):
        s1 = Stdout()
        self.assertTrue(True)

    def test_stdout_except(self):
        subtype = DataTypeWithSecondary()
        self.assertRaises(Exception, Stdout, subtype)


class TestParseTypes(unittest.TestCase):
    def test_parse_primitive_str(self):
        t = get_instantiated_type(str)
        self.assertIsInstance(t, String)

    def test_parse_type_str(self):
        t = get_instantiated_type(String)
        self.assertIsInstance(t, String)

    def test_parse_instantiated_str(self):
        t = get_instantiated_type(String())
        self.assertIsInstance(t, String)

    def test_parse_primitive_int(self):
        t = get_instantiated_type(int)
        self.assertIsInstance(t, Int)

    def test_parse_type_int(self):
        t = get_instantiated_type(Int)
        self.assertIsInstance(t, Int)

    def test_parse_instantiated_int(self):
        t = get_instantiated_type(Int())
        self.assertIsInstance(t, Int)

    def test_parse_primitive_float(self):
        t = get_instantiated_type(float)
        self.assertIsInstance(t, Float)

    def test_parse_type_float(self):
        t = get_instantiated_type(Float)
        self.assertIsInstance(t, Float)

    def test_parse_instantiated_float(self):
        t = get_instantiated_type(Float())
        self.assertIsInstance(t, Float)

    def test_parse_primitive_bool(self):
        t = get_instantiated_type(bool)
        self.assertIsInstance(t, Boolean)

    def test_parse_type_bool(self):
        t = get_instantiated_type(Boolean)
        self.assertIsInstance(t, Boolean)

    def test_parse_instantiated_bool(self):
        t = get_instantiated_type(Boolean())
        self.assertIsInstance(t, Boolean)
