import unittest
from typing import List, Dict, Any


from janis_core.tests.testtools import (
    SingleTestTool,
    ArrayTestTool,
    TestTool,
    TestToolWithSecondaryOutput,
    TestTypeWithSecondary,
    TestWorkflowWithStepInputExpression,
)

import cwlgen

import janis_core.translations.cwl as cwl
from janis_core import (
    WorkflowBuilder,
    ToolOutput,
    ToolInput,
    String,
    CommandTool,
    Stdout,
    InputSelector,
    Array,
    File,
    WildcardSelector,
    StringFormatter,
    ToolArgument,
)
from janis_core.tool.documentation import InputDocumentation
from janis_core.translations import CwlTranslator
from janis_core.types import CpuSelector, MemorySelector
from janis_core.workflow.workflow import InputNode


class TestCwlTypesConversion(unittest.TestCase):
    pass


class TestCwlMisc(unittest.TestCase):
    def test_str_tool(self):
        t = TestTool()
        self.assertEqual(cwl_testtool, t.translate("cwl", to_console=False))


class TestCwlTranslatorOverrides(unittest.TestCase):
    def setUp(self):
        self.translator = CwlTranslator()

    def test_stringify_WorkflowBuilder(self):
        cwlobj = cwlgen.Workflow("wid")
        expected = """\
#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.0
inputs: {}
outputs: {}
steps: {}
id: wid
"""
        self.assertEqual(
            expected, self.translator.stringify_translated_workflow(cwlobj)
        )

    def test_stringify_tool(self):
        cwlobj = cwlgen.CommandLineTool("tid")
        self.assertEqual(
            "#!/usr/bin/env cwl-runner\nclass: CommandLineTool\ncwlVersion: v1.0\nid: tid\n",
            self.translator.stringify_translated_tool(cwlobj),
        )

    def test_stringify_inputs(self):
        d = {"inp1": 1}
        self.assertEqual("inp1: 1\n", self.translator.stringify_translated_inputs(d))

    def test_workflow_filename(self):
        w = WorkflowBuilder("wid")
        self.assertEqual("wid.cwl", self.translator.workflow_filename(w))

    def test_tools_filename(self):
        self.assertEqual(
            "TestTranslationtool.cwl", self.translator.tool_filename(TestTool())
        )

    def test_inputs_filename(self):
        w = WorkflowBuilder("wid")
        self.assertEqual("wid-inp.yml", self.translator.inputs_filename(w))

    def test_resources_filename(self):
        w = WorkflowBuilder("wid")
        self.assertEqual("wid-resources.yml", self.translator.resources_filename(w))


class TestCwlArraySeparators(unittest.TestCase):
    # Based on https://www.commonwl.org/user_guide/09-array-inputs/index.html

    def test_regular_input_bindingin(self):
        t = ToolInput("filesA", Array(String()), prefix="-A", position=1)
        cwltoolinput = cwl.translate_tool_input(t)
        self.assertDictEqual(
            {
                "id": "filesA",
                "label": "filesA",
                "type": {"items": "string", "type": "array"},
                "inputBinding": {"prefix": "-A", "position": 1},
            },
            cwltoolinput.get_dict(),
        )

    def test_nested_input_binding(self):
        t = ToolInput(
            "filesB",
            Array(String()),
            prefix="-B=",
            separate_value_from_prefix=False,
            position=2,
            prefix_applies_to_all_elements=True,
        )
        cwltoolinput = cwl.translate_tool_input(t)
        self.assertDictEqual(
            {
                "id": "filesB",
                "label": "filesB",
                "type": {
                    "items": "string",
                    "type": "array",
                    "inputBinding": {"prefix": "-B=", "separate": False},
                },
                "inputBinding": {"position": 2},
            },
            cwltoolinput.get_dict(),
        )

    def test_separated_input_bindingin(self):
        t = ToolInput(
            "filesC",
            Array(String()),
            prefix="-C=",
            separate_value_from_prefix=False,
            position=4,
            separator=",",
        )
        cwltoolinput = cwl.translate_tool_input(t)
        self.assertDictEqual(
            {
                "id": "filesC",
                "label": "filesC",
                "type": {"items": "string", "type": "array"},
                "inputBinding": {
                    "prefix": "-C=",
                    "itemSeparator": ",",
                    "separate": False,
                    "position": 4,
                },
            },
            cwltoolinput.get_dict(),
        )

    def test_optional_array_prefixes(self):
        t = ToolInput(
            "filesD",
            Array(String(), optional=True),
            prefix="-D",
            prefix_applies_to_all_elements=True,
        )
        cwltoolinput = cwl.translate_tool_input(t)

        self.assertDictEqual(
            {
                "id": "filesD",
                "label": "filesD",
                "type": [
                    {
                        "inputBinding": {"prefix": "-D"},
                        "items": "string",
                        "type": "array",
                    },
                    "null",
                ],
            },
            cwltoolinput.get_dict(),
        )


class TestCwlSelectorsAndGenerators(unittest.TestCase):
    def test_input_selector_base(self):
        input_sel = InputSelector("random")
        self.assertEqual(
            "$(inputs.random)",
            cwl.translate_input_selector(input_sel, code_environment=False),
        )

    def test_input_selector_base_codeenv(self):
        input_sel = InputSelector("random")
        self.assertEqual(
            "inputs.random",
            cwl.translate_input_selector(input_sel, code_environment=True),
        )

    def test_input_value_none_codeenv(self):
        self.assertEqual(
            None, cwl.CwlTranslator.unwrap_expression(None, code_environment=True)
        )

    def test_input_value_none_nocodeenv(self):
        self.assertEqual(
            None, cwl.CwlTranslator.unwrap_expression(None, code_environment=False)
        )

    def test_input_value_string_codeenv(self):
        self.assertEqual(
            '"TestString"',
            cwl.CwlTranslator.unwrap_expression("TestString", code_environment=True),
        )

    def test_input_value_string_nocodeenv(self):
        self.assertEqual(
            "TestString",
            cwl.CwlTranslator.unwrap_expression("TestString", code_environment=False),
        )

    def test_input_value_int_codeenv(self):
        self.assertEqual(
            42, cwl.CwlTranslator.unwrap_expression(42, code_environment=True)
        )

    def test_input_value_int_nocodeenv(self):
        self.assertEqual(
            42, cwl.CwlTranslator.unwrap_expression(42, code_environment=False)
        )

    # def test_input_value_filename_codeenv(self):
    #     import uuid
    #     fn = Filename(guid=str(uuid.uuid4()))
    #     self.assertEqual(
    #         '"generated-" + Math.random().toString(16).substring(2, 8) + ""',
    #         cwl.get_input_value_from_potential_selector_or_generator(fn, code_environment=True)
    #     )
    #
    # def test_input_value_filename_nocodeenv(self):
    #     import uuid
    #     fn = Filename(guid=str(uuid.uuid4()))
    #     self.assertEqual(
    #         '$("generated-" + Math.random().toString(16).substring(2, 8) + "")',
    #         cwl.get_input_value_from_potential_selector_or_generator(fn, code_environment=False)
    #     )

    def test_input_value_inpselect_codeenv(self):
        inp = InputSelector("threads")
        self.assertEqual(
            "inputs.threads",
            cwl.CwlTranslator.unwrap_expression(inp, code_environment=True),
        )

    def test_input_value_inpselect_nocodeenv(self):
        inp = InputSelector("threads")
        self.assertEqual(
            "$(inputs.threads)",
            cwl.CwlTranslator.unwrap_expression(inp, code_environment=False),
        )

    def test_input_value_wildcard(self):
        self.assertRaises(
            Exception, cwl.CwlTranslator.unwrap_expression, value=WildcardSelector("*")
        )

    def test_input_value_cpuselect_codeenv(self):
        inp = CpuSelector()
        self.assertEqual(
            "inputs.runtime_cpu" "",
            cwl.CwlTranslator.unwrap_expression(inp, code_environment=True),
        )

    def test_input_value_cpuselect_nocodeenv(self):
        inp = CpuSelector()
        self.assertEqual(
            "$(inputs.runtime_cpu)",
            cwl.CwlTranslator.unwrap_expression(inp, code_environment=False),
        )

    def test_input_value_memselect_codeenv(self):
        inp = MemorySelector()
        self.assertEqual(
            "inputs.runtime_memory",
            cwl.CwlTranslator.unwrap_expression(inp, code_environment=True),
        )

    def test_input_value_memselect_nocodeenv(self):
        inp = MemorySelector()
        self.assertEqual(
            "$(inputs.runtime_memory)",
            cwl.CwlTranslator.unwrap_expression(inp, code_environment=False),
        )

    def test_input_value_cwl_callable(self):
        class NonCallableCwl:
            def cwl(self):
                return "unbelievable"

        self.assertEqual(
            "unbelievable", cwl.CwlTranslator.unwrap_expression(NonCallableCwl())
        )

    def test_input_value_cwl_noncallable(self):
        class NonCallableCwl:
            def __init__(self):
                self.cwl = None

        self.assertRaises(
            Exception,
            cwl.CwlTranslator.unwrap_expression,
            value=NonCallableCwl(),
            tool_id=None,
        )

    def test_string_formatter(self):
        b = StringFormatter("no format")
        res = cwl.CwlTranslator.unwrap_expression(b)
        self.assertEqual("no format", res)

    def test_string_formatter_one_string_param(self):
        b = StringFormatter("there's {one} arg", one="a string")
        res = cwl.CwlTranslator.unwrap_expression(b)
        self.assertEqual('$("there\'s {one} arg".replace(/\{one\}/g, "a string"))', res)

    def test_string_formatter_one_input_selector_param(self):
        b = StringFormatter("an input {arg}", arg=InputSelector("random_input"))
        res = cwl.CwlTranslator.unwrap_expression(b, code_environment=False)
        self.assertEqual(
            '$("an input {arg}".replace(/\{arg\}/g, inputs.random_input))', res
        )

    def test_string_formatter_two_param(self):
        # vardict input format
        b = StringFormatter(
            "{tumorName}:{normalName}",
            tumorName=InputSelector("tumorInputName"),
            normalName=InputSelector("normalInputName"),
        )
        res = cwl.CwlTranslator.unwrap_expression(b)
        self.assertEqual(
            '$("{tumorName}:{normalName}".replace(/\{tumorName\}/g, inputs.tumorInputName).replace(/\{normalName\}/g, inputs.normalInputName))',
            res,
        )

    def test_escaped_characters(self):
        trans = cwl.CwlTranslator
        translated = trans.translate_tool_internal(TestTool())
        arg: cwlgen.CommandLineBinding = translated.arguments[0]
        self.assertEqual('test:\\\\t:escaped:\\\\n:characters"', arg.valueFrom)


class TestCwlEnvVar(unittest.TestCase):
    def test_environment1(self):
        t = CwlTranslator().translate_tool_internal(tool=TestTool())
        envvar: cwlgen.EnvVarRequirement = [
            t for t in t.requirements if t._req_class == "EnvVarRequirement"
        ][0]
        envdef: cwlgen.EnvVarRequirement.EnvironmentDef = envvar.envDef[0]
        self.assertEqual("test1", envdef.envName)
        self.assertEqual("$(inputs.testtool)", envdef.envValue)


class TestCwlTranslateInput(unittest.TestCase):
    def test_translate_input(self):
        inp = InputNode(
            None,
            identifier="testIdentifier",
            datatype=String(),
            default="defaultValue",
            doc=InputDocumentation("docstring"),
            value=None,
        )
        tinp = cwl.translate_workflow_input(inp)

        self.assertEqual("testIdentifier", tinp.id)
        self.assertIsNone(tinp.label)
        self.assertIsNone(tinp.secondaryFiles)
        self.assertEqual("docstring", tinp.doc)
        self.assertIsNone(None, tinp.inputBinding)
        self.assertEqual("string", tinp.type)
        self.assertEqual("defaultValue", tinp.default)

    def test_secondary_file_translation(self):
        inp = InputNode(
            None,
            identifier="testIdentifier",
            datatype=TestTypeWithSecondary(),
            default=None,
            value=None,
        )
        tinp = cwl.translate_workflow_input(inp)

        self.assertEqual("File", tinp.type)
        self.assertListEqual(["^.txt"], tinp.secondaryFiles)


# PUT RIGHT HERE


class TestCwlGenerateInput(unittest.TestCase):
    def setUp(self):
        self.translator = cwl.CwlTranslator()

    def test_input_in_input_value_nooptional_nodefault(self):
        wf = WorkflowBuilder("test_cwl_input_in_input_value_nooptional_nodefault")
        wf.input("inpId", String(), default="1")
        self.assertDictEqual({"inpId": "1"}, self.translator.build_inputs_file(wf))

    def test_input_in_input_value_nooptional_default(self):
        wf = WorkflowBuilder("test_cwl_input_in_input_value_nooptional_default")
        wf.input("inpId", String(), default="1")
        self.assertDictEqual({"inpId": "1"}, self.translator.build_inputs_file(wf))

    def test_input_in_input_value_optional_nodefault(self):
        wf = WorkflowBuilder("test_cwl_input_in_input_value_optional_nodefault")
        wf.input("inpId", String(optional=True), default="1")
        self.assertDictEqual({"inpId": "1"}, self.translator.build_inputs_file(wf))

    def test_input_in_input_value_optional_default(self):
        wf = WorkflowBuilder("test_cwl_input_in_input_value_optional_default")
        wf.input("inpId", String(optional=True), default="1")
        self.assertDictEqual({"inpId": "1"}, self.translator.build_inputs_file(wf))

    def test_input_in_input_novalue_nooptional_nodefault(self):
        wf = WorkflowBuilder("test_cwl_input_in_input_novalue_nooptional_nodefault")
        wf.input("inpId", String())
        # included because no value, no default, and not optional
        # self.assertDictEqual({"inpId": None}, self.translator.build_inputs_file(wf))
        self.assertDictEqual({}, self.translator.build_inputs_file(wf))

    def test_input_in_input_novalue_nooptional_default(self):
        wf = WorkflowBuilder("test_cwl_input_in_input_novalue_nooptional_default")
        wf.input("inpId", String(), default="2")
        self.assertDictEqual({"inpId": "2"}, self.translator.build_inputs_file(wf))

    def test_input_in_input_novalue_optional_nodefault(self):
        wf = WorkflowBuilder("test_cwl_input_in_input_novalue_optional_nodefault")
        wf.input("inpId", String(optional=True))
        # self.assertDictEqual({'inpId': None}, self.translator.build_inputs_file(wf))
        self.assertDictEqual({}, self.translator.build_inputs_file(wf))

    def test_input_in_input_novalue_optional_default(self):
        wf = WorkflowBuilder("test_cwl_input_in_input_novalue_optional_default")
        wf.input("inpId", String(optional=True), default="2")
        self.assertDictEqual({"inpId": "2"}, self.translator.build_inputs_file(wf))


class TestCwlMaxResources(unittest.TestCase):
    def test_cores(self):
        tool = TestTool()
        resources = CwlTranslator.build_resources_input(tool.wrapped_in_wf(), {})
        self.assertEqual(2, resources["testtranslationtool_runtime_cpu"])

    def test_max_cores(self):
        tool = TestTool()
        resources = CwlTranslator.build_resources_input(
            tool.wrapped_in_wf(), {}, max_cores=1
        )
        self.assertEqual(1, resources["testtranslationtool_runtime_cpu"])

    def test_memory(self):
        tool = TestTool()
        resources = CwlTranslator.build_resources_input(tool.wrapped_in_wf(), {})
        self.assertEqual(2, resources["testtranslationtool_runtime_memory"])

    def test_max_memory(self):
        tool = TestTool()
        resources = CwlTranslator.build_resources_input(
            tool.wrapped_in_wf(), {}, max_mem=1
        )
        self.assertEqual(1, resources["testtranslationtool_runtime_memory"])


class TestEmptyContainer(unittest.TestCase):
    def test_empty_container_raises(self):

        self.assertRaises(
            Exception, CwlTranslator().translate_tool_internal, SingleTestTool()
        )

    def test_empty_container(self):
        c = CwlTranslator().translate_tool_internal(
            SingleTestTool(), allow_empty_container=True
        )
        self.assertNotIn("DockerRequirement", c.requirements)


class TestCwlSingleToMultipleInput(unittest.TestCase):
    def test_add_single_to_array_edge(self):
        w = WorkflowBuilder("test_add_single_to_array_edge")
        w.input("inp1", str)
        w.step("stp1", ArrayTestTool(inputs=w.inp1))

        c, _, _ = CwlTranslator().translate(
            w, to_console=False, allow_empty_container=True
        )
        self.assertEqual(cwl_multiinput, c)


class TestPackedWorkflow(unittest.TestCase):
    def test_simple(self):
        w = WorkflowBuilder("test_add_single_to_array_edge")
        w.step("ech", SingleTestTool(inputs="Hello"), doc="Print 'Hello'")
        c = CwlTranslator.translate_workflow_to_all_in_one(
            w, allow_empty_container=True
        )
        print(CwlTranslator.stringify_translated_workflow(c))


class TestCompleteOperators(unittest.TestCase):
    def test_step_input(self):

        ret, _, _ = TestWorkflowWithStepInputExpression().translate(
            "cwl", to_console=False
        )
        self.assertEqual(cwl_stepinput, ret)


cwl_testtool = """\
#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0
label: TestTranslationtool
requirements:
  DockerRequirement:
    dockerPull: ubuntu:latest
  EnvVarRequirement:
    envDef:
    - envName: test1
      envValue: $(inputs.testtool)
  InlineJavascriptRequirement: {}
  ShellCommandRequirement: {}
inputs:
- id: testtool
  label: testtool
  type: string
outputs:
- id: std
  label: std
  type: stdout
baseCommand: echo
arguments:
- position: 0
  valueFrom: test:\\\\t:escaped:\\\\n:characters"
id: TestTranslationtool
"""


cwl_multiinput = """\
#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.0
requirements:
  InlineJavascriptRequirement: {}
  MultipleInputFeatureRequirement: {}
  StepInputExpressionRequirement: {}
inputs:
  inp1:
    id: inp1
    type: string
outputs: {}
steps:
  stp1:
    in:
      inputs:
        id: inputs
        source:
        - inp1
        linkMerge: merge_nested
    run: tools/ArrayStepTool.cwl
    out:
    - outs
id: test_add_single_to_array_edge
"""

cwl_stepinput = """\
#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.0
label: 'TEST: WorkflowWithStepInputExpression'
requirements:
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}
inputs:
  mystring:
    id: mystring
    type:
    - string
    - 'null'
  mystring_backup:
    id: mystring_backup
    type:
    - string
    - 'null'
outputs:
  out:
    id: out
    type: File
    outputSource: print/out
steps:
  print:
    in:
      _print_inp_mystring:
        id: _print_inp_mystring
        source: mystring
      _print_inp_mystringbackup:
        id: _print_inp_mystringbackup
        source: mystring_backup
      inp:
        id: inp
        valueFrom: |-
          $((inputs._print_inp_mystring != null) ? inputs._print_inp_mystring : inputs._print_inp_mystringbackup)
    run: tools/EchoTestTool.cwl
    out:
    - out
id: TestWorkflowWithStepInputExpression
"""
