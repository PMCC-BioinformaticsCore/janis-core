import unittest
from typing import List, Optional

from janis_core.translations import CwlTranslator
from janis_core.types import String, Boolean, Float, Int, File, Array, Filename

from janis_core.code.pythontool import PythonTool

from janis_core.tool.tool import TOutput, TInput


class PythonEchoTool(PythonTool):
    @staticmethod
    def code_block(name: str, infile: Filename, flag: bool = True, testvalue="test"):
        """
        :param name: Name of the parameter
        :param infile: File to write to fout
        :param flag: Random boolean
        :param testvalue:
        :return:
        """
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

    def bind_metadata(self):
        return


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
        self.assertIsInstance(in2_infile.intype, Filename)
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

    def test_whole(self):
        out = PythonEchoTool().translate("wdl", to_console=False)
        print(out)
        # self.assertEqual(wdl, out)

    def test_whole2(self):
        test = CwlTranslator.translate_code_tool_internal(PythonEchoTool())
        print(test)

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


wdl = """\
version development

task echo_tool {
  input {
    Int? runtime_cpu
    Int? runtime_memory
    String name
    String? infile
    Boolean? flag
    String? testvalue
    String runtime_disks
  }
  command <<<
    
cat <<EOT >> echo_tool-script.py

import argparse, json, sys
from typing import Optional, List, Dict, Any
cli = argparse.ArgumentParser("Argument parser for Janis PythonTool")
cli.add_argument("--name", type=str, required=True)
cli.add_argument("--infile", type=str, help='File to write to fout')
cli.add_argument("--flag", action='store_true', help='Random boolean')
cli.add_argument("--testvalue", type=str)

String = str
Filename = str
Boolean = str
Int = int
Float = float
Double = float
File = str
Directory = str
Stdout = str
Stderr = str
Array = List
class PythonTool:
    File = str
    Directory = str



def code_block(name: str, infile: Filename, flag: bool = True, testvalue="test"):
    \"""
    :param name: Name of the parameter
    :param infile: File to write to fout
    :param flag: Random boolean
    :param testvalue:
    :return:
    \"""
    isset = "is set" if flag else str(type(flag))

    with open(infile, "w+") as tf:
        tf.write(f"Hello, {name}")

    return {"out": f"Hello, {name} ({isset})", "fout": infile}


try:
    args = cli.parse_args()
    result = code_block(name=args.name, infile=args.infile, flag=args.flag, testvalue=args.testvalue)
    print(json.dumps(result))
except Exception as e:
    print(str(e), file=sys.stderr)
    raise

EOT
    python echo_tool-script.py \\
      --name '~{name}' \\
      --infile '~{select_first([infile, "generated"])}' \\
      ~{if defined(select_first([flag, true])) then "--flag" else ""} \\
      ~{if defined(select_first([testvalue, "test"])) then ("--testvalue '" + select_first([testvalue, "test"]) + "'") else ""}
  >>>
  runtime {
    disks: runtime_disks
    docker: "python:3.8.1"
    memory: "~{select_first([runtime_memory, 4])}G"
    zones: "australia-southeast1-b"
  }
  output {
    String out = read_json(stdout())["out"]
    File fout = read_json(stdout())["fout"]
  }
}"""
