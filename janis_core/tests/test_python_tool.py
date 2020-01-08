import unittest
from typing import List, Optional

from janis_core.types import String, Boolean, Float, Int, File, Array

from janis_core.code.pythontool import PythonTool, get_instantiated_type_override

from janis_core.tool.tool import TOutput, TInput


class PythonEchoTool(PythonTool):
    @staticmethod
    def code_block(name: str, infile: Optional[PythonTool.File], flag2: bool):
        isset = "is set" if flag2 else str(type(flag2))

        with open(infile, "w+") as tf:
            tf.write(f"Hello, {name}")

        return {"out": f"Hello, {name} ({isset})", "fout": infile}

    def outputs(self) -> List[TOutput]:
        return [TOutput("out", String()), TOutput("fout", File())]

    def id(self) -> str:
        return "echo_tool"

    def version(self):
        return "v0.1.0"


# class PythonToolCodeBuilderTests(unittest.TestCase):
#     def test_inputs(self):
#         ins = PythonEchoTool().inputs()
#         print(len(ins))
#
#     def test_build_code_block(self):
#         script = PythonEchoTool().prepared_script()
#         print(script)
#
#     def test_generate_input_binding(self):
#         test1 = PythonTool.generate_cli_binding_for_input(TInput("test1", String()))
#         print(test1)
#
#     def test_generate_input_binding2(self):
#         test1 = PythonTool.generate_cli_binding_for_input(TInput("testflag", Boolean()))
#         print(test1)
#
#     def test_translation(self):
#         self.assertTrue(True)
#         # PythonEchoTool().translate("cwl")


class PythonToolTypeOverride(unittest.TestCase):
    def test_file(self):
        t = get_instantiated_type_override(PythonTool.File)
        self.assertIsInstance(t, File)
        self.assertFalse(t.optional)

    def test_optional_file(self):
        t = get_instantiated_type_override(Optional[PythonTool.File])
        self.assertIsInstance(t, File)
        self.assertTrue(t.optional)

    def test_list_file(self):
        t = get_instantiated_type_override(List[PythonTool.File])
        self.assertIsInstance(t, Array)
        self.assertIsInstance(t.subtype(), File)
        self.assertFalse(t.optional)
