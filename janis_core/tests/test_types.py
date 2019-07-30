import unittest

from janis_core import Array, String, Stdout, File


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
