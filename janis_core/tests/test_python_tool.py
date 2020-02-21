import unittest
from typing import List, Optional

from janis_core.types import String, Boolean, Float, Int, File, Array

from janis_core.code.pythontool import PythonTool

from janis_core.tool.tool import TOutput, TInput


class PythonEchoTool(PythonTool):
    @staticmethod
    def code_block(
        name: str, infile: Optional[File], flag: bool = True, testvalue="test"
    ):
        isset = "is set" if flag else str(type(flag))

        with open(infile, "w+") as tf:
            tf.write(f"Hello, {name}")

        return {"out": f"Hello, {name} ({isset})", "fout": infile}

    def outputs(self) -> List[TOutput]:
        return [TOutput("out", String()), TOutput("fout", File())]

    def id(self) -> str:
        return "echo_tool"

    def version(self):
        return "v0.1.0"


class PythonToolCodeBuilderTests(unittest.TestCase):
    def test_inputs(self):
        ins = PythonEchoTool().inputs()
        in1_name = ins[0]
        self.assertFalse(in1_name.intype.optional)
        self.assertEqual("name", in1_name.tag)
        self.assertIsInstance(in1_name.intype, String)
        self.assertIsNone(in1_name.default)

        in2_infile = ins[1]
        self.assertTrue(in2_infile.intype.optional)
        self.assertEqual("infile", in2_infile.tag)
        self.assertIsInstance(in2_infile.intype, File)
        self.assertIsNone(in2_infile.default)

        in3_flag = ins[2]
        self.assertTrue(in3_flag.intype.optional)
        self.assertEqual("flag", in3_flag.tag)
        self.assertIsInstance(in3_flag.intype, Boolean)
        self.assertTrue(in3_flag.default)

        in4_testvalue = ins[3]
        self.assertTrue(in4_testvalue.intype.optional)
        self.assertEqual("testvalue", in4_testvalue.tag)
        self.assertIsInstance(in4_testvalue.intype, String)
        self.assertEqual("test", in4_testvalue.default)

    # def test_build_code_block(self):
    #     script = PythonEchoTool().prepared_script()
    #     print(script)
    #
    # def test_generate_input_binding(self):
    #     test1 = PythonTool.generate_cli_binding_for_input(TInput("test1", String()))
    #     print(test1)
    #
    # def test_generate_input_binding2(self):
    #     test1 = PythonTool.generate_cli_binding_for_input(TInput("testflag", Boolean()))
    #     print(test1)
    #
    # def test_translation(self):
    #     self.assertTrue(True)
    #     PythonEchoTool().translate("cwl")
