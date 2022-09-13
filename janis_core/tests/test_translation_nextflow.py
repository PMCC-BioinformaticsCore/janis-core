import unittest
from typing import List, Dict, Any, Optional

from janis_core.operators.logical import If, IsDefined
from janis_core.operators.standard import ReadContents

from janis_core.tests.testtools import (
    SingleTestTool,
    ArrayTestTool,
    TestTool,
    TestToolV2,
    TestInputQualityTool,
    TestTypeWithSecondary,
    TestWorkflowWithStepInputExpression,
    EchoTestTool,
    FilenameGeneratedTool,
    OperatorResourcesTestTool,
    TestWorkflowWithConditionStep,
    TestWorkflowThatOutputsArraysOfSecondaryFiles,
    TestWorkflowWithAliasSelectorWorkflow,
    TestToolWithSecondaryInput,
    TestToolWithAppendedSecondaryOutput,
    TestToolWithReplacedSecondaryOutput,
    TestTypeWithAppendedSecondary,
    TestSplitTextTool,
    TestSumTool,
    TestJoinArrayTool,
    TestFileInput,
    TestFileWithSecondaryInput,
)


from janis_core.translations.nextflow import NextflowTranslator
from janis_core import (
    WorkflowBuilder,
    ToolInput,
    InputSelector,
    WildcardSelector,
    StringFormatter,
    CommandToolBuilder,
    ToolOutput,
    DataType,
    JoinOperator,
)
from janis_core.tool.documentation import InputDocumentation
from janis_core.translations import NextflowTranslator as translator
from janis_core.types import CpuSelector, MemorySelector, Stdout
from janis_core.workflow.workflow import InputNode

from janis_core.operators.standard import FirstOperator
from janis_core import Array, String, Stdout, File, Int, Float, Boolean, Filename


class DataTypeWithSecondary(File):
    @staticmethod
    def name() -> str:
        return "test_secondary"

    @staticmethod
    def secondary_files():
        return [".txt", ".csv"]


class DataTypeNoSecondary(File):
    @staticmethod
    def name() -> str:
        return "test_no_secondary"


class TestNextflowWfToolInputs(unittest.TestCase):
    def test_first_selector(self):

        workflow = TestWorkflowWithConditionStep()
        step_keys = list(workflow.step_nodes.keys())

        step_id = "print"
        tool = workflow.step_nodes[step_id].tool
        inputs = translator.gen_wf_tool_inputs(tool)
        expected = {"inp": "[$params.mystring, $get_string.out.out].first()"}

        self.assertEqual(expected, inputs)

    def test_simple(self):
        w1 = TestWorkflowThatOutputsArraysOfSecondaryFiles()
        w1_step_keys = list(w1.step_nodes.keys())

        expected = {"testtool": "$params.inp"}
        self.assertEqual(
            expected,
            translator.gen_wf_tool_inputs(w1.step_nodes["stp"].tool),
        )

    def test_with_expression(self):
        w2 = TestWorkflowWithStepInputExpression()
        w2_step_keys = list(w2.step_nodes.keys())

        expected = {
            "inp": "$params.mystring ? $params.mystring : $params.mystring_backup"
        }
        self.assertEqual(
            expected,
            translator.gen_wf_tool_inputs(w2.step_nodes["print"].tool),
        )

    def test_multi_steps(self):
        w3 = TestWorkflowWithAliasSelectorWorkflow()
        w3_step_keys = list(w3.step_nodes.keys())

        expected1 = {"testtool": "$params.inp"}
        self.assertEqual(
            expected1,
            translator.gen_wf_tool_inputs(w3.step_nodes["stp1"].tool),
        )

        expected2 = {"inp": "$stp1.out.out"}
        self.assertEqual(
            expected2,
            translator.gen_wf_tool_inputs(w3.step_nodes["stp2"].tool),
        )


class TestNextflowPrepareInputVars(unittest.TestCase):
    def test_secondary_files(self):
        inp = ToolInput("bam", DataTypeWithSecondary(), prefix="-I")

        res = translator.gen_input_var_definition(inp, "bam")
        expected = "apply_prefix(bam[0], '-I ', 'False')"
        self.assertEqual(res, expected)

    def test_array_with_secondary_files(self):

        inp = ToolInput("bams", Array(DataTypeWithSecondary()), prefix="-I")

        res = translator.gen_input_var_definition(inp, "bams")
        expected = "apply_prefix(get_primary_files(bams).join(' '), '-I ', 'False')"

        self.assertEqual(res, expected)


class TestGenerateWfToolOutputs(unittest.TestCase):
    w1 = TestWorkflowThatOutputsArraysOfSecondaryFiles()
    w2 = TestWorkflowWithStepInputExpression()
    w3 = TestWorkflowWithAliasSelectorWorkflow()

    def test_without_prefix(self):
        assert translator.gen_wf_tool_outputs(self.w1) == {"out": "stp.out.out"}
        assert translator.gen_wf_tool_outputs(self.w2) == {"out": "print.out.out"}
        assert translator.gen_wf_tool_outputs(self.w3) == {"out": "stp1.out.out"}

    def test_with_prefix(self):
        expected1 = {"out": "subworkflow_stp.out.out"}
        self.assertEqual(
            expected1, translator.gen_wf_tool_outputs(self.w1, "subworkflow_")
        )

        expected2 = {"out": "subworkflowprint.out.out"}
        self.assertEqual(
            expected2, translator.gen_wf_tool_outputs(self.w2, "subworkflow")
        )


class TestTranslateStringFormatter(unittest.TestCase):
    any_tool = TestTool()

    def test_string_formatter(self):
        b = StringFormatter("no format")
        res = translator.translate_string_formatter(b, self.any_tool)
        self.assertEqual("no format", res)

    def test_string_formatter_one_string_param(self):
        b = StringFormatter("there's {one} arg", one="a string")
        res = translator.translate_string_formatter(b, self.any_tool)
        self.assertEqual("there's ${'a string'} arg", res)

    def test_string_formatter_one_input_selector_param(self):
        b = StringFormatter("an input {arg}", arg=InputSelector("testtool"))
        res = translator.translate_string_formatter(
            b, self.any_tool, inputs_dict=self.any_tool.inputs_map()
        )
        self.assertEqual("an input ${testtool}", res)

    def test_string_formatter_two_param(self):
        tool = TestInputQualityTool()
        b = StringFormatter(
            "{username}:{password}",
            username=InputSelector("user"),
            password=InputSelector("static"),
        )
        res = translator.translate_string_formatter(
            b, tool, inputs_dict=tool.inputs_map()
        )
        self.assertEqual(
            "${user}:${static}",
            res,
        )

    def test_escaped_characters(self):
        tool = TestInputQualityTool()
        b = StringFormatter(
            "{username}\\t{password}",
            username=InputSelector("user"),
            password=InputSelector("static"),
        )
        res = translator.translate_string_formatter(
            b, tool, inputs_dict=tool.inputs_map()
        )
        self.assertEqual("${user}\\t${static}", res)

        res2 = translator.translate_string_formatter(
            b, tool, inputs_dict=tool.inputs_map(), in_shell_script=True
        )
        self.assertEqual("${user}\\\\t${static}", res2)

    def test_expression_arg(self):
        tool = TestTool()
        b = StringFormatter(
            "{name}:{items}",
            name=InputSelector("testtool"),
            items=JoinOperator(InputSelector("arrayInp"), separator=";"),
        )

        res = translator.translate_string_formatter(
            b, tool, inputs_dict=tool.inputs_map()
        )
        self.assertEqual("${testtool}:${arrayInp.join(';')}", res)


class TestNextflowGenerateInput(unittest.TestCase):
    def setUp(self):
        self.translator = translator

    def test_input_in_input_value_nooptional_nodefault(self):
        wf = WorkflowBuilder(
            "test_nf_input_in_input_value_nooptional_nodefault",
        )
        wf.input("inpId", String())
        ai = {"inpId": "2"}
        self.assertDictEqual(
            {"inpId": "2"}, self.translator.build_inputs_file(wf, additional_inputs=ai)
        )

    def test_input_in_input_value_nooptional_default(self):
        wf = WorkflowBuilder("test_nf_input_in_input_value_nooptional_default")
        wf.input("inpId", String(), default="1")
        ai = {"inpId": "2"}
        self.assertDictEqual(
            {"inpId": "2"}, self.translator.build_inputs_file(wf, additional_inputs=ai)
        )

    def test_input_in_input_value_optional_nodefault(self):
        wf = WorkflowBuilder("test_nf_input_in_input_value_optional_nodefault")
        wf.input("inpId", String(optional=True))
        ai = {"inpId": "2"}
        self.assertDictEqual(
            {"inpId": "2"}, self.translator.build_inputs_file(wf, additional_inputs=ai)
        )

    def test_input_in_input_value_optional_default(self):
        wf = WorkflowBuilder("test_nf_input_in_input_value_optional_default")
        wf.input("inpId", String(optional=True), default="1")
        ai = {"inpId": "2"}
        self.assertDictEqual(
            {"inpId": "2"}, self.translator.build_inputs_file(wf, additional_inputs=ai)
        )

    def test_input_in_input_novalue_nooptional_nodefault(self):
        wf = WorkflowBuilder("test_nf_input_in_input_novalue_nooptional_nodefault")
        wf.input("inpId", String())
        # included because no value, no default, and not optional
        self.assertDictEqual({"inpId": ""}, self.translator.build_inputs_file(wf))

    def test_input_in_input_novalue_nooptional_default(self):
        wf = WorkflowBuilder("test_nf_input_in_input_novalue_nooptional_default")
        wf.input("inpId", String(), default="1")
        self.assertDictEqual({"inpId": "1"}, self.translator.build_inputs_file(wf))

    def test_input_in_input_novalue_optional_nodefault(self):
        wf = WorkflowBuilder("test_nf_input_in_input_novalue_optional_nodefault")
        wf.input("inpId", String(optional=True))
        self.assertDictEqual({"inpId": ""}, self.translator.build_inputs_file(wf))

    def test_input_in_input_novalue_optional_default(self):
        wf = WorkflowBuilder("test_nf_input_in_input_novalue_optional_default")
        wf.input("inpId", String(optional=True), default="1")
        self.assertDictEqual({"inpId": "1"}, self.translator.build_inputs_file(wf))

    def test_bool_input_type_value(self):
        wf = WorkflowBuilder("test_bool_input")
        wf.input("inpBool", Boolean())
        ai = {"inpBool": True}
        self.assertDictEqual(
            ai, self.translator.build_inputs_file(wf, additional_inputs=ai)
        )

    def test_bool_input_type_novalue_default(self):
        wf = WorkflowBuilder("test_bool_input")
        wf.input("inpBool", Boolean(), default=False)
        self.assertDictEqual({"inpBool": False}, self.translator.build_inputs_file(wf))

    def test_bool_input_type_value_default(self):
        wf = WorkflowBuilder("test_bool_input")
        wf.input("inpBool", Boolean(), default=False)
        ai = {"inpBool": True}
        self.assertDictEqual(
            ai, self.translator.build_inputs_file(wf, additional_inputs=ai)
        )

    def test_file_type_with_secondary_value(self):
        wf = WorkflowBuilder("test_file_input")
        wf.input("inp", DataTypeWithSecondary())
        ai = {"inp": "/my/path/filename.doc"}
        self.assertDictEqual(
            {
                "inp": [
                    "/my/path/filename.doc",
                    "/my/path/filename.doc.txt",
                    "/my/path/filename.doc.csv",
                ]
            },
            self.translator.build_inputs_file(wf, additional_inputs=ai),
        )

    def test_file_type_no_secondary_value(self):
        wf = WorkflowBuilder("test_file_input")
        wf.input("inp", DataTypeNoSecondary())
        ai = {"inp": "/my/path/filename.doc"}
        self.assertDictEqual(
            {"inp": "/my/path/filename.doc"},
            self.translator.build_inputs_file(wf, additional_inputs=ai),
        )

    def test_file_type_with_secondary_no_value(self):
        wf = WorkflowBuilder("test_file_input")
        wf.input("inp", DataTypeWithSecondary())
        self.assertDictEqual(
            {"inp": f"/{translator.NO_FILE_PATH_PREFIX}1"},
            self.translator.build_inputs_file(wf),
        )

    def test_file_type_no_secondary_no_value(self):
        wf = WorkflowBuilder("test_file_input")
        wf.input("inp", DataTypeNoSecondary())
        self.assertDictEqual(
            {"inp": f"/{translator.NO_FILE_PATH_PREFIX}1"},
            self.translator.build_inputs_file(wf),
        )

    def test_file_type_no_secondary_multiple_inputs(self):
        wf = WorkflowBuilder("test_file_input")
        wf.input("inp", DataTypeNoSecondary())
        wf.input("inp2", DataTypeNoSecondary())
        self.assertDictEqual(
            {
                "inp": f"/{translator.NO_FILE_PATH_PREFIX}1",
                "inp2": f"/{translator.NO_FILE_PATH_PREFIX}2",
            },
            self.translator.build_inputs_file(wf),
        )

    def test_bool_filename_type_no_value(self):
        wf = WorkflowBuilder("test_filename_input")
        wf.input("inp", String(), default="somefancyname")
        wf.input(
            "inpFile",
            Filename(prefix=InputSelector("inp"), suffix="part1", extension=".doc"),
        )
        # Note: we don't want to set default filenames into json input file
        self.assertDictEqual(
            {"inp": "somefancyname", "inpFile": ""},
            self.translator.build_inputs_file(wf),
        )


class TestGenerateNfProcessForCommandTool(unittest.TestCase):
    maxDiff: Optional[int] = None
    def test_stdout_out_tool(self):
        p = translator.gen_process_from_cmdtool(EchoTestTool())
        p_str = p.get_string()
        expected = f"""
process EchoTestTool {{
  input:
    val inp

  output:
    path "${{'janisstdout_EchoTestTool'}}" , emit: out

  publishDir "$launchDir/EchoTestTool"
  memory "${{params.runtime_memory ? params.runtime_memory + 'GB': ''}}"
  cpus "${{params.runtime_cpu ? params.runtime_cpu : ''}}"
  disk "${{params.runtime_disks ? params.runtime_disks + 'GB': ''}}"
  time "${{params.runtime_seconds + 's'}}"

  script:
    def inpWithPrefix = apply_prefix(inp, ' ', 'False')
    def runtime_memory = params.runtime_memory
    def runtime_cpu = params.runtime_cpu
    def runtime_disks = params.runtime_disks
    def runtime_seconds = params.runtime_seconds
    
    \"\"\"
    echo \\
    $inpWithPrefix | tee janisstdout_EchoTestTool
    \"\"\"
}}


"""
        self.assertEqual(expected, p_str)

    def test_operator_resource_tool(self):
        p = translator.gen_process_from_cmdtool(OperatorResourcesTestTool())
        expected = f"""
process EchoTestTool
{{
  input:
    path inputFile
    val outputFiles

  output:
    path "${{'janisstdout_EchoTestTool'}}" , emit: out

  publishDir "$launchDir/EchoTestTool"
  memory "${{params.runtime_memory ? params.runtime_memory + 'GB': ''}}"
  cpus "${{params.runtime_cpu ? params.runtime_cpu : ''}}"
  disk "${{params.runtime_disks ? params.runtime_disks + 'GB': ''}}"
  time "${{params.runtime_seconds + 's'}}"

  script:

    def inputFileWithPrefix = apply_prefix(inputFile, ' ', 'False')

    def outputFilesWithPrefix = apply_prefix(outputFiles, ' ', 'False')

    def runtime_memory = params.runtime_memory

    def runtime_cpu = params.runtime_cpu

    def runtime_disks = params.runtime_disks

    def runtime_seconds = params.runtime_seconds
    \"\"\"
    echo \\
    $inputFileWithPrefix | tee janisstdout_EchoTestTool
    \"\"\"
}}


"""
        self.assertEqual(expected, p.get_string())

    def test_filename_generated_tool(self):
        p = translator.gen_process_from_cmdtool(FilenameGeneratedTool())

        expected = f"""
process filenamegeneratedtool
{{
  input:
    val inp
    val inpOptional
    path fileInp
    path fileInpOptional
    val generatedInp
    val generatedInpOptional
    val generatedFileInp
    val generatedFileInpOptional

  output:
    val "${{'*'}}" , emit: out

  publishDir "$launchDir/filenamegeneratedtool"
  memory "${{params.runtime_memory ? params.runtime_memory + 'GB': ''}}"
  cpus "${{params.runtime_cpu ? params.runtime_cpu : ''}}"
  disk "${{params.runtime_disks ? params.runtime_disks + 'GB': ''}}"
  time "${{params.runtime_seconds + 's'}}"

  script:

    def generatedInp = generatedInp && generatedInp != 'None' ? generatedInp : inp + '' + ''

    def generatedInpOptional = generatedInpOptional && generatedInpOptional != 'None' ? generatedInpOptional : inpOptional + '' + ''

    def generatedFileInp = generatedFileInp && generatedFileInp != 'None' ? generatedFileInp : fileInp.simpleName + '.transformed' + '.fnp'

    def generatedFileInpOptional = generatedFileInpOptional && generatedFileInpOptional != 'None' ? generatedFileInpOptional : fileInpOptional.simpleName + '.optional' + '.txt'

    def inpWithPrefix = apply_prefix(inp, ' ', 'False')

    def inpOptionalWithPrefix = optional(inpOptional, ' ', 'False')

    def fileInpWithPrefix = apply_prefix(fileInp, ' ', 'False')

    def fileInpOptionalWithPrefix = optional(fileInpOptional, ' ', 'False')

    def generatedInpWithPrefix = optional(generatedInp, ' ', 'False')

    def generatedInpOptionalWithPrefix = optional(generatedInpOptional, ' ', 'False')

    def generatedFileInpWithPrefix = optional(generatedFileInp, ' ', 'False')

    def generatedFileInpOptionalWithPrefix = optional(generatedFileInpOptional, ' ', 'False')

    def runtime_memory = params.runtime_memory

    def runtime_cpu = params.runtime_cpu

    def runtime_disks = params.runtime_disks

    def runtime_seconds = params.runtime_seconds
    \"\"\"
    echo \\
    $generatedInpWithPrefix \\
    $generatedInpOptionalWithPrefix \\
    $generatedFileInpWithPrefix \\
    $generatedFileInpOptionalWithPrefix | tee janisstdout_filenamegeneratedtool
    \"\"\"
}}


"""
        self.assertEqual(expected, p.get_string())

    def test_tool_with_secondary_input(self):
        p = translator.gen_process_from_cmdtool(
            TestToolWithSecondaryInput()
        )

        expected = f"""
process CatTestTool
{{
  input:
    path inp

  output:
    path "${{'janisstdout_CatTestTool'}}" , emit: out

  publishDir "$launchDir/CatTestTool"
  memory "${{params.runtime_memory ? params.runtime_memory + 'GB': ''}}"
  cpus "${{params.runtime_cpu ? params.runtime_cpu : ''}}"
  disk "${{params.runtime_disks ? params.runtime_disks + 'GB': ''}}"
  time "${{params.runtime_seconds + 's'}}"

  script:

    def inpWithPrefix = apply_prefix(inp[0], ' ', 'False')

    def runtime_memory = params.runtime_memory

    def runtime_cpu = params.runtime_cpu

    def runtime_disks = params.runtime_disks

    def runtime_seconds = params.runtime_seconds
    \"\"\"
    cat \\
    $inpWithPrefix | tee janisstdout_CatTestTool
    \"\"\"
}}


"""
        self.assertEqual(expected, p.get_string())

    def test_tool_with_secondary_output(self):
        p = translator.gen_process_from_cmdtool(
            TestToolWithAppendedSecondaryOutput()
        )

        expected = f"""
process TestTranslationtool
{{
  input:
    val testtool
    val arrayInp

  output:
    tuple path("${{testtool + '.bam'}}"), path("${{testtool + '.bam.bai'}}") , emit: out

  publishDir "$launchDir/TestTranslationtool"
  memory "${{params.runtime_memory ? params.runtime_memory + 'GB': ''}}"
  cpus "${{params.runtime_cpu ? params.runtime_cpu : ''}}"
  disk "${{params.runtime_disks ? params.runtime_disks + 'GB': ''}}"
  time "${{params.runtime_seconds + 's'}}"

  script:

    def testtoolWithPrefix = apply_prefix(testtool, ' ', 'False')

    def arrayInpWithPrefix = optional(arrayInp.join(' '), ' ', 'False')

    def runtime_memory = params.runtime_memory

    def runtime_cpu = params.runtime_cpu

    def runtime_disks = params.runtime_disks

    def runtime_seconds = params.runtime_seconds
    \"\"\"
    echo \\
    test:\\t:escaped:\\n:characters" | tee janisstdout_TestTranslationtool
    \"\"\"
}}


"""
        self.assertEqual(expected, p.get_string())

    def test_tool_with_replaced_secondary_output(self):
        p = translator.gen_process_from_cmdtool(
            TestToolWithReplacedSecondaryOutput()
        )
        expected = f"""
process TestTranslationtool
{{
  input:
    val testtool
    val arrayInp

  output:
    tuple path("${{testtool + '.bam'}}"), path("${{testtool + '.bai'}}") , emit: out

  publishDir "$launchDir/TestTranslationtool"
  memory "${{params.runtime_memory ? params.runtime_memory + 'GB': ''}}"
  cpus "${{params.runtime_cpu ? params.runtime_cpu : ''}}"
  disk "${{params.runtime_disks ? params.runtime_disks + 'GB': ''}}"
  time "${{params.runtime_seconds + 's'}}"

  script:

    def testtoolWithPrefix = apply_prefix(testtool, ' ', 'False')

    def arrayInpWithPrefix = optional(arrayInp.join(' '), ' ', 'False')

    def runtime_memory = params.runtime_memory

    def runtime_cpu = params.runtime_cpu

    def runtime_disks = params.runtime_disks

    def runtime_seconds = params.runtime_seconds
    \"\"\"
    echo \\
    test:\\t:escaped:\\n:characters" | tee janisstdout_TestTranslationtool
    \"\"\"
}}


"""

        self.assertEqual(expected, p.get_string())


class TestGenerateNfProcessForPythonCodeTool(unittest.TestCase):
    def test_str_input(self):
        p = translator.gen_process_from_codetool(TestSplitTextTool())
        expected = f"""
process TestSplitTextTool
{{
  input:
    path PYTHON_CODE_FILE_PATH
    val inp

  output:
    val "${{file("$workDir/janis_out_out").text.replace('[', '').replace(']', '').split(', ')}}" , emit: out

  publishDir "$launchDir/TestSplitTextTool"
  memory "${{params.runtime_memory ? params.runtime_memory + 'GB': ''}}"
  cpus "${{params.runtime_cpu ? params.runtime_cpu : ''}}"
  disk "${{params.runtime_disks ? params.runtime_disks + 'GB': ''}}"
  time "${{params.runtime_seconds + 's'}}"

  script:
    \"\"\"
    #!/usr/bin/env python
    from TestSplitTextTool import code_block
    import os
    import json

    result = code_block(inp="$inp")

    work_dir = os.getenv("PYENV_DIR")
    for key in result:
        with open(os.path.join("$workDir", f"janis_out_{{key}}"), "w") as f:
            f.write(json.dumps(result[key]))
    \"\"\"
}}


"""
        self.assertEqual(expected, p.get_string())

    def test_int_input(self):
        p = translator.gen_process_from_codetool(TestSumTool())
        expected = f"""
process TestSumTool
{{
  input:
    path PYTHON_CODE_FILE_PATH
    val inp1
    val inp2

  output:
    val "${{file("$workDir/janis_out_out").text.replace('[', '').replace(']', '').split(', ')}}" , emit: out

  publishDir "$launchDir/TestSumTool"
  memory "${{params.runtime_memory ? params.runtime_memory + 'GB': ''}}"
  cpus "${{params.runtime_cpu ? params.runtime_cpu : ''}}"
  disk "${{params.runtime_disks ? params.runtime_disks + 'GB': ''}}"
  time "${{params.runtime_seconds + 's'}}"

  script:
    \"\"\"
    #!/usr/bin/env python
    from TestSumTool import code_block
    import os
    import json

    result = code_block(inp1=$inp1, inp2=$inp2)

    work_dir = os.getenv("PYENV_DIR")
    for key in result:
        with open(os.path.join("$workDir", f"janis_out_{{key}}"), "w") as f:
            f.write(json.dumps(result[key]))
    \"\"\"
}}


"""
        self.assertEqual(expected, p.get_string())

    def test_array_input(self):
        p = translator.gen_process_from_codetool(TestJoinArrayTool())
        expected = f"""
process TestJoinArrayTool
{{
  input:
    path PYTHON_CODE_FILE_PATH
    val inp

  output:
    val "${{file("$workDir/janis_out_out").text.replace('[', '').replace(']', '').split(', ')}}" , emit: out

  publishDir "$launchDir/TestJoinArrayTool"
  memory "${{params.runtime_memory ? params.runtime_memory + 'GB': ''}}"
  cpus "${{params.runtime_cpu ? params.runtime_cpu : ''}}"
  disk "${{params.runtime_disks ? params.runtime_disks + 'GB': ''}}"
  time "${{params.runtime_seconds + 's'}}"

  script:
    \"\"\"
    #!/usr/bin/env python
    from TestJoinArrayTool import code_block
    import os
    import json

    result = code_block(inp="$inp".split(" "))

    work_dir = os.getenv("PYENV_DIR")
    for key in result:
        with open(os.path.join("$workDir", f"janis_out_{{key}}"), "w") as f:
            f.write(json.dumps(result[key]))
    \"\"\"
}}


"""
        self.assertEqual(expected, p.get_string())

    def test_file_input(self):
        p = translator.gen_process_from_codetool(TestFileInput())
        expected = f"""
process TestFileInput
{{
  input:
    path PYTHON_CODE_FILE_PATH
    path inp

  output:
    val "${{file("$workDir/janis_out_out").text.replace('[', '').replace(']', '').split(', ')}}" , emit: out

  publishDir "$launchDir/TestFileInput"
  memory "${{params.runtime_memory ? params.runtime_memory + 'GB': ''}}"
  cpus "${{params.runtime_cpu ? params.runtime_cpu : ''}}"
  disk "${{params.runtime_disks ? params.runtime_disks + 'GB': ''}}"
  time "${{params.runtime_seconds + 's'}}"

  script:
    \"\"\"
    #!/usr/bin/env python
    from TestFileInput import code_block
    import os
    import json

    result = code_block(inp="$inp")

    work_dir = os.getenv("PYENV_DIR")
    for key in result:
        with open(os.path.join("$workDir", f"janis_out_{{key}}"), "w") as f:
            f.write(json.dumps(result[key]))
    \"\"\"
}}


"""
        self.assertEqual(expected, p.get_string())

    def test_file_with_secondary_input(self):
        p = translator.gen_process_from_codetool(
            TestFileWithSecondaryInput()
        )
        expected = f"""
process TestFileWithSecondaryInput
{{
  input:
    path PYTHON_CODE_FILE_PATH
    path inp

  output:
    val "${{file("$workDir/janis_out_out").text.replace('[', '').replace(']', '').split(', ')}}" , emit: out

  publishDir "$launchDir/TestFileWithSecondaryInput"
  memory "${{params.runtime_memory ? params.runtime_memory + 'GB': ''}}"
  cpus "${{params.runtime_cpu ? params.runtime_cpu : ''}}"
  disk "${{params.runtime_disks ? params.runtime_disks + 'GB': ''}}"
  time "${{params.runtime_seconds + 's'}}"

  script:
    \"\"\"
    #!/usr/bin/env python
    from TestFileWithSecondaryInput import code_block
    import os
    import json

    result = code_block(inp="$inp".split(" ")[0])

    work_dir = os.getenv("PYENV_DIR")
    for key in result:
        with open(os.path.join("$workDir", f"janis_out_{{key}}"), "w") as f:
            f.write(json.dumps(result[key]))
    \"\"\"
}}


"""
        self.assertEqual(expected, p.get_string())
