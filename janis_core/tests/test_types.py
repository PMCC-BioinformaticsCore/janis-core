import unittest
from typing import Optional, List, Union

from janis_core import Array, String, Stdout, File, Int, Float, Boolean
from janis_core.types import get_instantiated_type, get_from_python_type


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
        self.assertEqual({"type": "array", "items": "string"}, d.save())

    def test_array_of_array_of_strings(self):
        ar = Array(Array(String()))
        d = ar.cwl_type()
        self.assertDictEqual(
            {"type": "array", "items": {"type": "array", "items": "string"}}, d.save()
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

    def test_parse_optional_file(self):
        t = get_instantiated_type(Optional[File])
        self.assertTrue(t.optional)
        self.assertIsInstance(t, File)

    def test_parse_union_type(self):
        t = get_instantiated_type(Union[str, int])
        self.assertIsInstance(t.subtypes[0], String)
        self.assertIsInstance(t.subtypes[1], Int)
        self.assertEqual("Union<String, Integer>", t.id())

    def test_parse_union_optional_types(self):
        t = get_instantiated_type(Union[Optional[str], int])
        self.assertTrue(t.optional)
        self.assertIsInstance(t.subtypes[0], String)
        self.assertIsInstance(t.subtypes[1], Int)
        self.assertEqual("Union<String, Integer>", t.id())


class TestPythonAnnotations(unittest.TestCase):
    def test_string(self):
        t = get_from_python_type(str)
        self.assertIsInstance(t, String)

    def test_1(self):
        t = get_from_python_type(Optional[bool])
        self.assertTrue(t.optional)
        self.assertIsInstance(t, Boolean)

    def test_2(self):
        t = get_from_python_type(List[str])
        self.assertFalse(t.optional)
        self.assertTrue(t.is_array())
        st = t.subtype()
        self.assertFalse(st.optional)
        self.assertIsInstance(st, String)

    def test_3(self):
        t = get_from_python_type(List[List[str]])
        self.assertFalse(t.optional)
        self.assertTrue(t.is_array())
        st = t.subtype()
        self.assertFalse(st.optional)
        self.assertTrue(st.is_array())
        sst = st.subtype()
        self.assertFalse(sst.optional)
        self.assertIsInstance(sst, String)

    def test_4(self):
        # Might be a fun python 3.6 thing here...
        t = get_from_python_type(List[Optional[str]])
        self.assertFalse(t.optional)
        self.assertTrue(t.is_array())
        st = t.subtype()
        self.assertTrue(st.optional)
        self.assertIsInstance(st, String)
