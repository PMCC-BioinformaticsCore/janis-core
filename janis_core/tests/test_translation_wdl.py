import unittest
from tempfile import TemporaryDirectory
from typing import Optional

from janis_core.deps import wdlgen
from janis_core.types import UnionType

import janis_core.translations.wdl as wdl
from janis_core import (
    WorkflowBuilder,
    ToolOutput,
    ToolInput,
    ToolArgument,
    String,
    CommandTool,
    Stdout,
    InputSelector,
    Array,
    File,
    Filename,
    WildcardSelector,
    Boolean,
    Int,
    Float,
    CommandToolBuilder,
)
from janis_core.operators import CpuSelector, StringFormatter
from janis_core.operators.logical import If, IsDefined
from janis_core.operators.standard import JoinOperator, ReadContents
from janis_core.tests.testtools import (
    TestTypeWithSecondary,
    TestTypeWithNonEscapedSecondary,
    EchoTestTool,
    SingleTestTool,
    FilenameGeneratedTool,
    TestTool,
    TestToolV2,
    TestToolWithSecondaryInput,
    TestWorkflowWithStepInputExpression,
    ArrayTestTool,
    OperatorResourcesTestTool,
    TestWorkflowThatOutputsArraysOfSecondaryFiles,
    TestForEach,
)
from janis_core.translations import WdlTranslator
from janis_core.utils.scatter import ScatterDescription, ScatterMethod


class MultipleEcho(CommandTool):
    @staticmethod
    def tool():
        return "TestStepTool"

    @staticmethod
    def base_command():
        return "echo"

    def inputs(self):
        return [
            ToolInput("input1", TxtSecondary, position=0),
            ToolInput("input2", String(optional=True), position=1),
            ToolInput("input3", String(optional=True), position=2),
            ToolInput("input4", String(optional=True), position=3),
        ]

    def friendly_name(self):
        return None

    def outputs(self):
        return [ToolOutput("out", Stdout)]

    @staticmethod
    def container():
        return None

    @staticmethod
    def version():
        return None


class TxtSecondary(File):
    def __init__(self, optional=False):
        super().__init__(optional=optional, extension=".txt")

    @staticmethod
    def secondary_files():
        return [".qt"]


class TestToolWithSecondaryOutput(TestTool):
    def outputs(self):
        return [
            ToolOutput(
                "out",
                TestTypeWithNonEscapedSecondary(),
                selector=InputSelector("testtool") + "/out",
            )
        ]


class TestTypeWithAlternateAndSecondary(File):
    def __init__(self, optional=False):
        super().__init__(
            optional=optional, extension=".txt", alternate_extensions={".text"}
        )

    @staticmethod
    def secondary_files():
        return ["^.file"]


class TestWdl(unittest.TestCase):
    def test_optional_array(self):
        t = Array(File(), optional=True)
        wdl = t.wdl()
        self.assertIsInstance(wdl, wdlgen.WdlType)
        self.assertTrue(wdl.optional)
        self.assertEqual("Array[File]?", wdl.get_string())


class TestWdlTranslatorOverrides(unittest.TestCase):
    def setUp(self):
        self.translator = WdlTranslator()

    def test_stringify_workflow(self):
        wdlobj = wdlgen.Workflow("wid", version="development")
        self.assertEqual(
            "version development\n\n\n\nworkflow wid {\n\n}",
            self.translator.stringify_translated_workflow(wdlobj),
        )

    def test_stringify_tool(self):
        wdlobj = wdlgen.Task("tid", version="development")
        self.assertEqual(
            "version development\n\ntask tid {\n\n}",
            self.translator.stringify_translated_tool(wdlobj),
        )

    def test_stringify_inputs(self):
        d = {"wid.inp1": 1}
        self.assertEqual(
            '{\n    "wid.inp1": 1\n}', self.translator.stringify_translated_inputs(d)
        )

    def test_workflow_filename(self):
        w = WorkflowBuilder("wid")
        self.assertEqual("wid.wdl", self.translator.workflow_filename(w))

    def test_tools_filename(self):
        self.assertEqual(
            "TestTranslationtool.wdl", self.translator.tool_filename(TestTool().id())
        )

    def test_inputs_filename(self):
        w = WorkflowBuilder("wid")
        self.assertEqual("wid-inp.json", self.translator.inputs_filename(w))

    def test_resources_filename(self):
        w = WorkflowBuilder("wid")
        self.assertEqual("wid-resources.json", self.translator.resources_filename(w))


class TestWdlTranslatorBuilders(unittest.TestCase):
    def test_inputs_generator_secondary_files(self):
        w = WorkflowBuilder("tst")
        w.input("wsec", TestTypeWithSecondary, default="test.ext")
        # w._add_input(Input("wsec", TestTypeWithSecondary(), value="test.ext"))
        inpsdict = WdlTranslator().build_inputs_file(w, merge_resources=False)
        self.assertEqual("test.ext", inpsdict.get("tst.wsec"))
        self.assertEqual("test.txt", inpsdict.get("tst.wsec_txt"))

    def test_inputs_generator_array_of_secondary_files(self):
        w = WorkflowBuilder("tst")
        w.input("wsec", Array(TestTypeWithSecondary()), default=["test.ext"])
        inpsdict = WdlTranslator().build_inputs_file(w, merge_resources=False)
        self.assertListEqual(["test.ext"], inpsdict.get("tst.wsec"))
        self.assertListEqual(["test.txt"], inpsdict.get("tst.wsec_txt"))

    def test_translate_single_to_array_edge(self):
        w = WorkflowBuilder("wf")
        w.input("inp", str)
        stp1 = w.step("stp1", TestTool(testtool=w.inp), ignore_missing=True)
        stp2 = w.step(
            "stp2", TestTool(arrayInp=stp1.std, testtool=w.inp), ignore_missing=True
        )

        outp = wdl.translate_step_node(stp2, stp2.id(), {}, set(), None)
        self.assertEqual(
            f"arrayInp=[{stp1.id()}.std]", outp.get_string().split("\n")[3].strip()
        )


class TestWdlSelectorsAndGenerators(unittest.TestCase):
    def test_input_selector_base_stringenv(self):
        ti = {"random": ToolInput("random", String())}
        input_sel = InputSelector("random")
        self.assertEqual(
            "~{random}",
            wdl.translate_input_selector(input_sel, ti, string_environment=True),
        )

    def test_input_selector_base_nostringenv(self):
        ti = {"random": ToolInput("random", String())}
        input_sel = InputSelector("random")
        self.assertEqual(
            "random",
            wdl.translate_input_selector(input_sel, ti, string_environment=False),
        )

    def test_input_value_none_stringenv(self):
        self.assertEqual(
            "", wdl.WdlTranslator.unwrap_expression(None, None, string_environment=True)
        )

    def test_input_value_none_nostringenv(self):
        self.assertEqual(
            "",
            wdl.WdlTranslator.unwrap_expression(None, None, string_environment=False),
        )

    def test_input_value_string_stringenv(self):
        self.assertEqual(
            "TestString",
            wdl.WdlTranslator.unwrap_expression(
                "TestString", None, string_environment=True
            ),
        )

    def test_input_value_string_nostringenv(self):
        self.assertEqual(
            '"TestString"',
            wdl.WdlTranslator.unwrap_expression(
                "TestString", None, string_environment=False
            ),
        )

    def test_input_value_int_stringenv(self):
        self.assertEqual(
            str(42),
            wdl.WdlTranslator.unwrap_expression(42, None, string_environment=True),
        )

    def test_input_value_int_nostringenv(self):
        self.assertEqual(
            str(42),
            wdl.WdlTranslator.unwrap_expression(42, None, string_environment=False),
        )

    def test_input_value_filename_stringenv(self):

        fn = Filename()
        self.assertEqual(
            fn.generated_filename(),
            wdl.WdlTranslator.unwrap_expression(fn, None, string_environment=True),
        )

    def test_input_value_filename_nostringenv(self):

        fn = Filename()
        self.assertEqual(
            '"%s"' % fn.generated_filename(),
            wdl.WdlTranslator.unwrap_expression(fn, None, string_environment=False),
        )

    def test_input_value_wildcard(self):
        self.assertRaises(
            Exception,
            wdl.WdlTranslator.unwrap_expression,
            value=WildcardSelector("*"),
            tool_id=None,
        )

    # def test_input_value_cpuselect_stringenv(self):
    #     # CpuSelector relies on their being a runtime_cpu attribute,
    #     # this test will assume it's present, and '' will test ensure
    #     # that it is actually present
    #     ti = {"runtime_cpu": ToolInput("runtime_cpu", Int(), default=1)}
    #     inp = CpuSelector()
    #     self.assertEqual(
    #         "~{select_first([runtime_cpu, 1])}",
    #         wdl.WdlTranslator.unwrap_expression(inp, ti, string_environment=True),
    #     )

    # def test_input_value_cpuselect_nostringenv(self):
    #     # CpuSelector relies on their being a runtime_cpu attribute,
    #     # this test will assume it's present, and '' will test ensure
    #     # that it is actually present
    #
    #     ti = {"runtime_cpu": ToolInput("runtime_cpu", Int(), default=1)}
    #     inp = CpuSelector()
    #     self.assertEqual(
    #         "select_first([runtime_cpu, 1])",
    #         wdl.WdlTranslator.unwrap_expression(inp, ti, string_environment=False),
    #     )

    def test_tool_input_value_default_cpuselect(self):
        ti = ToolInput("threads", Int(), default=CpuSelector(), prefix="-t")
        tid = {"threads": ti}

        tr = wdl.translate_command_input(ti)
        self.assertEqual(
            "-t ~{select_first([threads, select_first([runtime_cpu, 1])])}",
            tr.get_string(),
        )

    def test_tool_input_value_default_cpuselect_nodefault(self):
        ti = ToolInput("threads", Int(), default=CpuSelector(None), prefix="-t")
        tid = {"threads": ti}

        tr = wdl.translate_command_input(ti)
        self.assertEqual("-t ~{select_first([threads, runtime_cpu])}", tr.get_string())

    def test_tool_input_optional_array(self):
        ti = ToolInput(
            "adapter",
            input_type=Array(String(), optional=True),
            prefix="-a",
            prefix_applies_to_all_elements=True,
        )
        tr = wdl.translate_command_input(ti)
        self.assertEqual(
            '~{if (defined(adapter) && length(select_first([adapter])) > 0) then "-a \'" + sep("\' -a \'", select_first([adapter])) + "\'" else ""}',
            tr.get_string(),
        )

    # def test_input_value_memselect_stringenv(self):
    #     inp = MemorySelector()
    #     self.assertEqual(
    #         "~{floor(runtime_memory)}",
    #         wdl.WdlTranslator.unwrap_expression(inp, string_environment=True)
    #     )
    #
    # def test_input_value_memselect_nostringenv(self):
    #     inp = MemorySelector()
    #     self.assertEqual(
    #         "floor(runtime_memory)",
    #         wdl.WdlTranslator.unwrap_expression(inp, string_environment=False)
    #     )

    def test_input_value_wdl_callable(self):
        class CallableWdl:
            def wdl(self):
                return "unbelievable"

        self.assertEqual(
            "unbelievable", wdl.WdlTranslator.unwrap_expression(CallableWdl(), None)
        )

    def test_input_value_wdl_noncallable(self):
        class NonCallableWdl:
            def __init__(self):
                self.wdl = None

        self.assertRaises(
            Exception,
            wdl.WdlTranslator.unwrap_expression,
            value=NonCallableWdl(),
            tool_id=None,
        )

    def test_string_formatter(self):
        b = StringFormatter("no format")
        res = wdl.WdlTranslator.unwrap_expression(b, None, string_environment=True)
        self.assertEqual("no format", res)

    def test_string_formatter_one_string_param(self):
        b = StringFormatter("there's {one} arg", one="a string")
        res = wdl.WdlTranslator.unwrap_expression(b, None, string_environment=True)
        self.assertEqual("there's a string arg", res)

    def test_string_formatter_one_input_selector_param(self):
        d = {"random_input": ToolInput("random_input", String())}
        b = StringFormatter("an input {arg}", arg=InputSelector("random_input"))
        res = wdl.WdlTranslator.unwrap_expression(b, d, string_environment=True)
        self.assertEqual("an input ~{random_input}", res)

    def test_string_formatter_two_param(self):
        # vardict input format
        d = {
            "tumorInputName": ToolInput("tumorInputName", String()),
            "normalInputName": ToolInput("normalInputName", String()),
        }
        b = StringFormatter(
            "{tumorName}:{normalName}",
            tumorName=InputSelector("tumorInputName"),
            normalName=InputSelector("normalInputName"),
        )
        res = wdl.WdlTranslator.unwrap_expression(b, d, string_environment=True)
        self.assertEqual("~{tumorInputName}:~{normalInputName}", res)

    def test_escaped_characters(self):
        trans = wdl.WdlTranslator
        translated = trans.translate_tool_internal(TestTool())
        arg = translated.command[-1].arguments[0]
        self.assertEqual("'test:\\t:escaped:\\n:characters\"'", arg.value)

    # test removed as optional placeholders don't throw errors anymore
    # def test_string_formatter_optional_inpselect_no_default(self):
    #     # will throw
    #     ti = {"ti": ToolInput("ti", String(optional=True))}
    #     b = StringFormatter("{place} michael", place=InputSelector("ti"))
    #     self.assertRaises(Exception, wdl.WdlTranslator.unwrap_expression, b, ti)

    def test_string_formatter_optional_inpselect_with_default(self):
        ti = {"ti": ToolInput("ti", String(optional=True), default="hi")}
        b = StringFormatter("{place} michael", place=InputSelector("ti"))
        res = wdl.WdlTranslator.unwrap_expression(b, ti, string_environment=True)
        self.assertEqual('~{select_first([ti, "hi"])} michael', res)

    def test_resolve_filename_in_inpselect(self):
        fn = Filename(extension=".ext")
        ti = {"ti": ToolInput("ti", fn)}
        b = StringFormatter("fn: {place}", place=InputSelector("ti"))
        res = wdl.WdlTranslator.unwrap_expression(b, ti)
        self.assertEqual(
            f'"fn: ~{{select_first([ti, "{fn.generated_filename()}"])}}"', res
        )


class TestWDLFilenameGeneration(unittest.TestCase):
    def test_1(self):
        tool = FilenameGeneratedTool()
        mapped = [
            a.get_string()
            for a in WdlTranslator.build_command_from_inputs(tool.inputs())
        ]

        self.assertEqual("'~{select_first([generatedInp, \"~{inp}\"])}'", mapped[0])
        self.assertEqual(
            "'~{select_first([generatedInpOptional, \"~{inpOptional}\"])}'", mapped[1]
        )
        self.assertEqual(
            '\'~{select_first([generatedFileInp, "~{basename(fileInp, ".txt")}.transformed.fnp"])}\'',
            mapped[2],
        )
        self.assertEqual(
            '\'~{select_first([generatedFileInpOptional, "~{basename(fileInpOptional, ".txt")}.optional.txt"])}\'',
            mapped[3],
        )


class TestWdlGenerateInput(unittest.TestCase):
    def setUp(self):
        self.translator = wdl.WdlTranslator()

    def test_input_in_input_value_nooptional_nodefault(self):
        wf = WorkflowBuilder("test_input_in_inputfile")
        wf.input("inpId", String(), default="1")

        self.assertDictEqual(
            {"test_input_in_inputfile.inpId": "1"},
            self.translator.build_inputs_file(wf),
        )

    def test_input_in_input_value_nooptional_default(self):
        wf = WorkflowBuilder("test_input_in_inputfile")
        wf.input("inpId", String(), default="1")

        self.assertDictEqual(
            {"test_input_in_inputfile.inpId": "1"},
            self.translator.build_inputs_file(wf),
        )

    def test_input_in_input_value_optional_nodefault(self):
        wf = WorkflowBuilder("test_input_in_inputfile")
        wf.input("inpId", String(optional=True), default="1")

        self.assertDictEqual(
            {"test_input_in_inputfile.inpId": "1"},
            self.translator.build_inputs_file(wf),
        )

    def test_input_in_input_value_optional_default(self):
        wf = WorkflowBuilder("test_input_in_inputfile")
        wf.input("inpId", String(optional=True), default="1")

        self.assertDictEqual(
            {"test_input_in_inputfile.inpId": "1"},
            self.translator.build_inputs_file(wf),
        )

    def test_input_in_input_novalue_nooptional_nodefault(self):
        wf = WorkflowBuilder("test_input_in_inputfile")
        wf.input("inpId", String())

        # included because no value, no default, and not optional
        self.assertDictEqual(
            {"test_input_in_inputfile.inpId": None},
            self.translator.build_inputs_file(wf),
        )

    def test_input_in_input_novalue_nooptional_default(self):
        wf = WorkflowBuilder("test_input_in_inputfile")
        wf.input("inpId", String(), default="2")

        # new interpretation: defaults appear in inputs
        self.assertDictEqual(
            {"test_input_in_inputfile.inpId": "2"},
            self.translator.build_inputs_file(wf),
        )
        # self.assertDictEqual({}, self.translator.build_inputs_file(wf))

    def test_overrided_input_optional_nodefault(self):
        wf = WorkflowBuilder("test_input_in_inputfile")
        wf.input("inpId", String(optional=True))

        ad = {"inpId": "2"}

        # new interpretation: defaults appear in inputs
        self.assertDictEqual(
            {"test_input_in_inputfile.inpId": "2"},
            self.translator.build_inputs_file(wf, additional_inputs=ad),
        )

    def test_overrided_input_optional_default(self):
        wf = WorkflowBuilder("test_input_in_inputfile")
        wf.input("inpId", String(optional=True), default="2")

        ad = {"inpId": "4"}

        # new interpretation: defaults appear in inputs
        self.assertDictEqual(
            {"test_input_in_inputfile.inpId": "4"},
            self.translator.build_inputs_file(wf, additional_inputs=ad),
        )

    def test_input_in_input_novalue_optional_nodefault(self):
        wf = WorkflowBuilder("test_input_in_inputfile")
        wf.input("inpId", String(optional=True))

        self.assertDictEqual({}, self.translator.build_inputs_file(wf))

    def test_input_in_input_novalue_optional_default(self):
        wf = WorkflowBuilder("test_input_in_inputfile")
        wf.input("inpId", String(optional=True), default="2")

        # new interpretation: defaults appear in inputs
        self.assertDictEqual(
            {"test_input_in_inputfile.inpId": "2"},
            self.translator.build_inputs_file(wf),
        )
        # self.assertDictEqual({}, self.translator.build_inputs_file(wf))

    def test_tool_output_with_input_selector(self):

        tool = TestToolWithSecondaryOutput()
        toolout = tool.outputs()[0]
        inmap = {t.id(): t for t in tool.inputs()}
        os = WdlTranslator.translate_tool_outputs([toolout], inmap, tool=tool)

        self.assertEqual("out", os[0].name)
        self.assertEqual('(testtool + "/out")', os[0].expression)

        self.assertEqual("out_txt", os[1].name)
        self.assertEqual('(testtool + "/out") + ".txt"', os[1].expression)

    def test_optional_tool_output_with_secondary(self):
        tool = TestToolWithSecondaryOutput()
        toolout = ToolOutput(
            "out",
            TestTypeWithNonEscapedSecondary(optional=True),
            selector=InputSelector("testtool"),
        )
        inmap = {t.id(): t for t in tool.inputs()}
        os = WdlTranslator.translate_tool_outputs([toolout], inmap, tool=tool)
        self.assertEqual(
            'File? out_txt = if defined(testtool) then (testtool + ".txt") else None',
            os[1].get_string(),
        )


class TestWdlToolInputGeneration(unittest.TestCase):
    def test_nodefault_nooptional_position(self):
        ti = ToolInput("tag", String(), position=0)
        resp = wdl.translate_command_input(ti)
        self.assertEqual("'~{tag}'", resp.get_string())

    def test_nodefault_nooptional_prefix_sep(self):
        ti = ToolInput("tag", String(), prefix="--amazing")
        resp = wdl.translate_command_input(ti)
        self.assertEqual("--amazing '~{tag}'", resp.get_string())

    def test_nodefault_nooptional_prefix_nosep(self):
        ti = ToolInput(
            "tag", String(), prefix="--amazing=", separate_value_from_prefix=False
        )
        resp = wdl.translate_command_input(ti)
        self.assertEqual("--amazing='~{tag}'", resp.get_string())

    def test_nodefault_optional_position(self):
        ti = ToolInput("tag", String(optional=True), position=0)
        resp = wdl.translate_command_input(ti)
        self.assertEqual(
            '~{if defined(tag) then ("\'" + tag + "\'") else ""}', resp.get_string()
        )

    def test_nodefault_optional_prefix_sep(self):
        ti = ToolInput("tag", String(optional=True), prefix="--amazing")
        resp = wdl.translate_command_input(ti)
        self.assertEqual(
            '~{if defined(tag) then ("--amazing \'" + tag + "\'") else ""}',
            resp.get_string(),
        )

    def test_nodefault_optional_prefix_nosep(self):
        ti = ToolInput(
            "tag",
            String(optional=True),
            prefix="--amazing=",
            separate_value_from_prefix=False,
        )
        resp = wdl.translate_command_input(ti)
        self.assertEqual(
            '~{if defined(tag) then ("--amazing=\'" + tag + "\'") else ""}',
            resp.get_string(),
        )

    def test_default_nooptional_position(self):
        # this will get turned into an optional
        ti = ToolInput("tag", String(), position=0, default="defval")
        resp = wdl.translate_command_input(ti)
        self.assertEqual("'~{select_first([tag, \"defval\"])}'", resp.get_string())

    def test_default_nooptional_prefix_sep(self):
        ti = ToolInput("tag", String(), prefix="--amazing", default="defval")
        resp = wdl.translate_command_input(ti)
        self.assertEqual(
            "--amazing '~{select_first([tag, \"defval\"])}'", resp.get_string()
        )

    def test_default_nooptional_prefix_nosep(self):
        ti = ToolInput(
            "tag",
            String(),
            prefix="--amazing=",
            separate_value_from_prefix=False,
            default="defval",
        )
        resp = wdl.translate_command_input(ti)
        self.assertEqual(
            "--amazing='~{select_first([tag, \"defval\"])}'", resp.get_string()
        )

    def test_default_optional_position(self):
        ti = ToolInput("tag", String(optional=True), position=0, default="defval")
        resp = wdl.translate_command_input(ti)
        self.assertEqual(
            '~{if defined(select_first([tag, "defval"])) then ("\'" + select_first([tag, "defval"]) + "\'") else ""}',
            resp.get_string(),
        )

    def test_default_optional_prefix_sep(self):
        ti = ToolInput(
            "tag", String(optional=True), prefix="--amazing", default="defval"
        )
        resp = wdl.translate_command_input(ti)
        self.assertEqual(
            '~{if defined(select_first([tag, "defval"])) then ("--amazing \'" + select_first([tag, "defval"]) + "\'") else ""}',
            resp.get_string(),
        )

    def test_default_optional_prefix_nosep(self):
        ti = ToolInput(
            "tag",
            String(optional=True),
            prefix="--amazing=",
            separate_value_from_prefix=False,
            default="defval",
        )
        resp = wdl.translate_command_input(ti)
        self.assertEqual(
            '~{if defined(select_first([tag, "defval"])) then ("--amazing=\'" + select_first([tag, "defval"]) + "\'") else ""}',
            resp.get_string(),
        )

    def test_bind_boolean(self):
        ti = ToolInput("tag", Boolean, prefix="--amazing", default=True)
        resp = wdl.translate_command_input(ti).get_string()
        self.assertEqual('~{if tag then "--amazing" else ""}', resp)

    def test_bind_optional_oolean_as_default(self):
        ti = ToolInput("tag", Boolean(optional=True), prefix="--amazing", default=True)
        resp = wdl.translate_command_input(ti).get_string()
        self.assertEqual(
            '~{if select_first([tag, true]) then "--amazing" else ""}', resp
        )

    def test_bind_boolean(self):
        ti = ToolInput("tag", Boolean(optional=True), prefix="--amazing")
        resp = wdl.translate_command_input(ti).get_string()
        self.assertEqual(
            '~{if (defined(tag) && select_first([tag])) then "--amazing" else ""}', resp
        )

    def test_array_prefix_each_element_non_quoted(self):
        ti = ToolInput(
            "tag", Array(Int), prefix="-i", prefix_applies_to_all_elements=True
        )
        resp = wdl.translate_command_input(ti).get_string()
        self.assertEqual(
            '~{if length(tag) > 0 then sep(" ", prefix("-i ", tag)) else ""}', resp
        )


class TestWdlInputTranslation(unittest.TestCase):
    def test_string_nooptional_nodefault(self):
        s = String()
        self.assertEqual("String", s.wdl(has_default=False).get_string())

    def test_string_nooptional_default(self):
        s = String()
        # As of 2019-07-10, the defaults are applied within the command input, so these can be null
        self.assertEqual("String?", s.wdl(has_default=True).get_string())

    def test_string_optional_nodefault(self):
        s = String(optional=True)
        self.assertEqual("String?", s.wdl(has_default=False).get_string())

    def test_string_optional_default(self):
        s = String(optional=True)
        self.assertEqual("String?", s.wdl(has_default=True).get_string())


class TestWdlEnvVar(unittest.TestCase):
    def test_environment1(self):
        t = WdlTranslator().translate_tool_internal(tool=TestTool())
        s = t.get_string()
        print(s)


class TestWdlMaxResources(unittest.TestCase):
    def test_cores(self):
        tool = TestTool()
        resources = WdlTranslator.build_resources_input(
            tool.wrapped_in_wf(), {}, is_root=True
        )
        self.assertEqual(
            2, resources["TestTranslationtoolWf.testtranslationtool_runtime_cpu"]
        )

    def test_max_cores(self):
        tool = TestTool()
        resources = WdlTranslator.build_resources_input(
            tool.wrapped_in_wf(), {}, max_cores=1, is_root=True
        )
        self.assertEqual(
            1, resources["TestTranslationtoolWf.testtranslationtool_runtime_cpu"]
        )

    def test_memory(self):
        tool = TestTool()
        resources = WdlTranslator.build_resources_input(
            tool.wrapped_in_wf(), {}, is_root=True
        )
        self.assertEqual(
            2, resources["TestTranslationtoolWf.testtranslationtool_runtime_memory"]
        )

    def test_max_memory(self):
        tool = TestTool()
        resources = WdlTranslator.build_resources_input(
            tool.wrapped_in_wf(), {}, max_mem=1, is_root=True
        )
        self.assertEqual(
            1, resources["TestTranslationtoolWf.testtranslationtool_runtime_memory"]
        )


class TestWdlScatterByMultipleFields(unittest.TestCase):
    def test_scatter_single(self):
        w = WorkflowBuilder("sbmf")
        w.input("inp", Array(str))
        w.input("inp2", str)

        step = w.step(
            "dotTool",
            SingleTestTool(input1=w.inp, input2=w.inp2),
            scatter=ScatterDescription(fields=["input1"], method=ScatterMethod.dot),
        )

        outp = wdl.translate_step_node(
            step, "A.SingleTestTool", {}, {"inp", "inp2"}, None
        )
        expected = """\
scatter (i in inp) {
   call A.SingleTestTool as dotTool {
    input:
      input1=i,
      input2=inp2
  }
}"""
        self.assertEqual(expected, outp.get_string(indent=0))

    def test_scatter_single_no_description(self):
        w = WorkflowBuilder("sbmf")
        w.input("inp", Array(str))
        w.input("inp2", str)

        step = w.step(
            "dotTool", SingleTestTool(input1=w.inp, input2=w.inp2), scatter="input1"
        )

        outp = wdl.translate_step_node(
            step, "A.SingleTestTool", {}, {"inp", "inp2"}, None
        )
        expected = """\
scatter (i in inp) {
   call A.SingleTestTool as dotTool {
    input:
      input1=i,
      input2=inp2
  }
}"""
        self.assertEqual(expected, outp.get_string(indent=0))

    def test_dot_2(self):
        w = WorkflowBuilder("sbmf")
        w.input("inp", Array(str))
        w.input("inp2", Array(str))

        step = w.step(
            "dotTool",
            SingleTestTool(input1=w.inp, input2=w.inp2),
            scatter=ScatterDescription(
                fields=["input1", "input2"], method=ScatterMethod.dot
            ),
        )

        outp = wdl.translate_step_node(
            step, "A.SingleTestTool", {}, {"inp", "inp2"}, None
        )
        expected = """\
scatter (Q in zip(inp, inp2)) {
   call A.SingleTestTool as dotTool {
    input:
      input1=Q.left,
      input2=Q.right
  }
}"""
        self.assertEqual(expected, outp.get_string(indent=0))

    def test_dot_3(self):
        w = WorkflowBuilder("sbmf")
        w.input("inp", Array(str))
        w.input("inp2", Array(str))
        w.input("inp3", Array(str))

        step = w.step(
            "dotTool",
            SingleTestTool(input1=w.inp, input2=w.inp2, input3=w.inp3),
            scatter=ScatterDescription(
                fields=["input1", "input2", "input3"], method=ScatterMethod.dot
            ),
        )

        outp = wdl.translate_step_node(
            step, "A.SingleTestTool", {}, {"inp", "inp2", "inp3"}, None
        )
        expected = """\
scatter (Q in zip(inp, zip(inp2, inp3))) {
   call A.SingleTestTool as dotTool {
    input:
      input1=Q.left,
      input2=Q.right.left,
      input3=Q.right.right
  }
}"""
        self.assertEqual(expected, outp.get_string(indent=0))

    def test_dot_4(self):
        w = WorkflowBuilder("sbmf")
        w.input("inp", Array(str))
        w.input("inp2", Array(str))
        w.input("inp3", Array(str))
        w.input("inp4", Array(str))

        step = w.step(
            "dotTool",
            SingleTestTool(input1=w.inp, input2=w.inp2, input3=w.inp3, input4=w.inp4),
            scatter=ScatterDescription(
                fields=["input1", "input2", "input3", "input4"],
                method=ScatterMethod.dot,
            ),
        )

        outp = wdl.translate_step_node(
            step, "A.SingleTestTool", {}, {"inp", "inp2", "inp3", "inp4"}, None
        )
        expected = """\
scatter (Q in zip(inp, zip(inp2, zip(inp3, inp4)))) {
   call A.SingleTestTool as dotTool {
    input:
      input1=Q.left,
      input2=Q.right.left,
      input3=Q.right.right.left,
      input4=Q.right.right.right
  }
}"""
        self.assertEqual(expected, outp.get_string(indent=0))

    def test_dot_2_secondary(self):
        w = WorkflowBuilder("sbmf")
        w.input("inp", Array(TxtSecondary))
        w.input("inp2", Array(str))

        step = w.step(
            "dotTool",
            MultipleEcho(input1=w.inp, input2=w.inp2),
            scatter=ScatterDescription(
                fields=["input1", "input2"], method=ScatterMethod.dot
            ),
        )

        outp = wdl.translate_step_node(
            step, "A.SingleTestTool", {}, {"inp", "inp2"}, None
        )
        expected = """\
scatter (Q in zip(transpose([inp, inp_qt]), inp2)) {
   call A.SingleTestTool as dotTool {
    input:
      input1=Q.left[0],
      input1_qt=Q.left[1],
      input2=Q.right
  }
}"""
        self.assertEqual(expected, outp.get_string(indent=0))


class TestRuntimeOverrideGenerator(unittest.TestCase):
    def test_basic(self):
        w = WorkflowBuilder("wb")
        w.input("inp", str)
        w.step("echo", SingleTestTool(input1=w.inp))
        w.step("echo_2", SingleTestTool(input1=w.inp))

        wf, _, _ = w.translate(
            "wdl",
            to_console=False,
            with_resource_overrides=True,
            allow_empty_container=True,
        )
        _tooldef = """\
workflow wb {
  input {
    String inp
    Int? echo_runtime_memory
    Int? echo_runtime_cpu
    Int? echo_runtime_disks
    Int? echo_runtime_seconds
    Int? echo_2_runtime_memory
    Int? echo_2_runtime_cpu
    Int? echo_2_runtime_disks
    Int? echo_2_runtime_seconds
  }
  call T.TestStepTool as echo {
    input:
      input1=inp,
      runtime_memory=echo_runtime_memory,
      runtime_cpu=echo_runtime_cpu,
      runtime_disks=echo_runtime_disks,
      runtime_seconds=echo_runtime_seconds
  }
  call T.TestStepTool as echo_2 {
    input:
      input1=inp,
      runtime_memory=echo_2_runtime_memory,
      runtime_cpu=echo_2_runtime_cpu,
      runtime_disks=echo_2_runtime_disks,
      runtime_seconds=echo_2_runtime_seconds
  }
}"""
        self.assertEqual(_tooldef, "\n".join(wf.split("\n")[4:]))


class TestLinkStatements(unittest.TestCase):
    def test_1(self):
        import janis_core as j

        class FileWithSec(j.File):
            def __init__(self, optional=False):
                super().__init__(optional=optional, extension=".txt")

            def secondary_files(self):
                return [".sec"]

        Tool = j.CommandToolBuilder(
            tool="ls",
            base_command=["ls"],
            inputs=[
                j.ToolInput(
                    "inp", FileWithSec, secondaries_present_as={".sec": "^.sec"}
                )
            ],
            outputs=[
                j.ToolOutput("std", j.Stdout),
                j.ToolOutput(
                    "out",
                    FileWithSec,
                    secondaries_present_as={".sec": "^.sec"},
                    glob=j.InputSelector("inp"),
                ),
            ],
            container="ubuntu:latest",
            version="v0.1.0",
        )

        Tool.translate("wdl")


class WorkflowWdlInputDefaultOperator(unittest.TestCase):
    def test_string_formatter(self):
        wf = WorkflowBuilder("wf")
        wf.input("sampleName", str)
        wf.input("platform", str)

        wf.input(
            "readGroupHeaderLine",
            String(),
            default=StringFormatter(
                "@RG\\tID:{name}\\tSM:{name}\\tLB:{name}\\tPL:{pl}",
                name=InputSelector("sampleName"),
                pl=InputSelector("platform"),
            ),
        )
        wf.step("print", EchoTestTool(inp=wf.readGroupHeaderLine))
        wf.output("out", source=wf.print)
        derived, _, _ = wf.translate("wdl", to_console=False)
        expected = """\
version development

import "tools/EchoTestTool_TEST.wdl" as E

workflow wf {
  input {
    String sampleName
    String platform
    String? readGroupHeaderLine = "@RG\\tID:~{sampleName}\\tSM:~{sampleName}\\tLB:~{sampleName}\\tPL:~{platform}"
  }
  call E.EchoTestTool as print {
    input:
      inp=select_first([readGroupHeaderLine, "@RG\\tID:~{sampleName}\\tSM:~{sampleName}\\tLB:~{sampleName}\\tPL:~{platform}"])
  }
  output {
    File out = print.out
  }
}"""
        self.assertEqual(expected, derived)


class TestWdlContainerOverride(unittest.TestCase):
    def test_tool_dict_override(self):
        expected_container = "container/override"

        tool = SingleTestTool()
        translated = tool.translate(
            "wdl", to_console=False, container_override={tool.id(): expected_container}
        )

        line = translated.splitlines()[-9].strip()
        self.assertEqual(f'docker: "{expected_container}"', line)

    def test_tool_string_override(self):
        expected_container = "container/override"

        tool = SingleTestTool()
        translated = tool.translate(
            "wdl", to_console=False, container_override=expected_container
        )

        line = translated.splitlines()[-9].strip()
        self.assertEqual(f'docker: "{expected_container}"', line)

    def test_tool_override_casecheck(self):
        expected_container = "container/override"

        tool = SingleTestTool()

        # Assert that our tool id is not UPPER, so when we override with the
        toolid_upper = tool.id().upper()
        self.assertNotEqual(tool.id(), toolid_upper)
        translated = tool.translate(
            "wdl",
            to_console=False,
            container_override={toolid_upper: expected_container},
        )

        line = translated.splitlines()[-9].strip()
        self.assertEqual(f'docker: "{expected_container}"', line)


class TestWDLRunRefs(unittest.TestCase):
    def test_two_similar_tools(self):
        w = WorkflowBuilder("testTwoToolsWithSameId")

        w.input("inp", str)
        w.step("stp1", TestTool(testtool=w.inp))
        w.step("stp2", TestToolV2(testtool=w.inp))

        wf_wdl, _ = WdlTranslator.translate_workflow(w)

        expected = """\
version development

import "tools/TestTranslationtool.wdl" as T
import "tools/TestTranslationtool_v0_0_2.wdl" as T2

workflow testTwoToolsWithSameId {
  input {
    String inp
  }
  call T.TestTranslationtool as stp1 {
    input:
      testtool=inp
  }
  call T2.TestTranslationtool as stp2 {
    input:
      testtool=inp
  }
}"""

        self.assertEqual(expected, wf_wdl.get_string())


class TestWdlSecondaryTranslation(unittest.TestCase):
    def test_secondary_connection(self):
        wf = WorkflowBuilder("wf")
        wf.input("ref", TestTypeWithSecondary)

        wf.step("stp", TestToolWithSecondaryInput(inp=wf.ref))

        wdlwf, _, _ = wf.translate("wdl", to_console=False)

        expected = """\
version development

import "tools/CatTestTool_TEST.wdl" as C

workflow wf {
  input {
    File ref
    File ref_txt
  }
  call C.CatTestTool as stp {
    input:
      inp=ref,
      inp_txt=ref_txt
  }
}"""
        self.assertEqual(expected, wdlwf)

    def test_array_secondary_connection(self):
        wf = WorkflowBuilder("wf")
        wf.input("ref", Array(TestTypeWithSecondary))

        wf.step("stp", TestToolWithSecondaryInput(inp=wf.ref), scatter="inp")

        wdlwf, _, _ = wf.translate("wdl", to_console=False)

        expected = """\
version development

import "tools/CatTestTool_TEST.wdl" as C

workflow wf {
  input {
    Array[File] ref
    Array[File] ref_txt
  }
  scatter (r in transpose([ref, ref_txt])) {
     call C.CatTestTool as stp {
      input:
        inp=r[0],
        inp_txt=r[1]
    }
  }
}"""
        self.assertEqual(expected, wdlwf)

    def test_workflow_secondary_outputs(self):
        wf = TestWorkflowThatOutputsArraysOfSecondaryFiles()
        wfwdl, _ = WdlTranslator.translate_workflow(wf)

        outs = [o.get_string() for o in wfwdl.outputs]
        self.assertEqual("Array[File] out = stp.out", outs[0])
        self.assertEqual("Array[File] out_txt = stp.out_txt", outs[1])

    def test_tool_with_secondary_and_alternates(self):
        tool = CommandToolBuilder(
            tool="test_secondary_and_alternates",
            base_command="cat",
            inputs=[
                ToolInput(
                    "inp",
                    TestTypeWithAlternateAndSecondary(),
                    position=1,
                    localise_file=True,
                )
            ],
            outputs=[
                ToolOutput(
                    "out",
                    TestTypeWithAlternateAndSecondary(),
                    selector=InputSelector("inp"),
                )
            ],
            container="ubtunu",
            version="TEST",
        )

        out = tool.translate("wdl", to_console=False)
        lines = out.splitlines(keepends=False)[-4:-2]
        l1 = "File out = basename(inp)"
        l2 = 'File out_file = sub(sub(basename(inp), "\\\\.txt$", ".file"), "\\\\.text$", ".file")'

        self.assertEqual(l1, lines[0].strip())
        self.assertEqual(l2, lines[1].strip())


class TestWDLCreateFilesAndDirectories(unittest.TestCase):

    initial_params = {
        "tool": "testCreateFilesAndDirectries",
        "version": "DEV",
        "container": "ubuntu",
        "base_command": "cat",
        "inputs": [ToolInput("inp", File), ToolInput("name", str)],
        "outputs": [ToolOutput("out", Stdout)],
    }

    def test_create_single_directory(self):
        command = CommandToolBuilder(
            **self.initial_params, directories_to_create="test-directory"
        )
        commands = WdlTranslator.build_commands_for_file_to_create(command)
        self.assertEqual(1, len(commands))

        self.assertEqual("mkdir -p 'test-directory'", commands[0].command)

    def test_create_single_directory_from_selector(self):
        command = CommandToolBuilder(
            **self.initial_params, directories_to_create=InputSelector("name")
        )
        commands = WdlTranslator.build_commands_for_file_to_create(command)
        self.assertEqual(1, len(commands))
        self.assertEqual("mkdir -p '~{name}'", commands[0].command)

    def test_create_single_directory_from_operator(self):
        command = CommandToolBuilder(
            **self.initial_params, directories_to_create=InputSelector("name") + "-out"
        )
        commands = WdlTranslator.build_commands_for_file_to_create(command)
        self.assertEqual(1, len(commands))
        self.assertEqual("mkdir -p '~{(name + \"-out\")}'", commands[0].command)

    def test_create_single_file_from_operator(self):
        command = CommandToolBuilder(
            **self.initial_params,
            files_to_create=[("my-path.txt", InputSelector("inp").contents())],
        )
        commands = WdlTranslator.build_commands_for_file_to_create(command)
        self.assertEqual(1, len(commands))
        expected = """\
cat <<EOT >> 'my-path.txt'
~{read_string(inp)}
EOT"""
        self.assertEqual(expected, commands[0].command)

    def test_create_single_file_path_from_operator(self):
        command = CommandToolBuilder(
            **self.initial_params,
            files_to_create=[
                (
                    StringFormatter("{name}.txt", name=InputSelector("name")),
                    "this is contents",
                )
            ],
        )
        commands = WdlTranslator.build_commands_for_file_to_create(command)
        self.assertEqual(1, len(commands))
        expected = """\
cat <<EOT >> '~{name}.txt'
this is contents
EOT"""
        self.assertEqual(expected, commands[0].command)


class TestCompleteOperators(unittest.TestCase):
    def test_list_operators(self):
        exp = WdlTranslator.unwrap_expression([1, 2, "three"])
        self.assertEqual('[1, 2, "three"]', exp)

    def test_step_input(self):

        ret, _, _ = TestWorkflowWithStepInputExpression().translate(
            "wdl", to_console=False
        )
        expected = """\
version development

import "tools/EchoTestTool_TEST.wdl" as E

workflow TestWorkflowWithStepInputExpression {
  input {
    String? mystring
    String? mystring_backup
  }
  call E.EchoTestTool as print {
    input:
      inp=if (defined(mystring)) then mystring else mystring_backup
  }
  output {
    File out = print.out
  }
}"""
        self.assertEqual(expected, ret)

    # def test_separator(self):
    #     tf = CommandToolBuilder(
    #         tool="test_sep_operator",
    #         base_command="echo",
    #         inputs=[ToolInput("inp", Array(String))],
    #         arguments=[
    #             ToolArgument(JoinOperator(InputSelector("inp"), ","), position=0)
    #         ],
    #         outputs=[ToolOutput("out", Stdout)],
    #         container="ubuntu:latest",
    #         version="v",
    #     )
    #
    #     tf.translate("wdl", to_disk=True, export_path="~/Desktop/tmp/wdltests/")

    def test_array_step_input(self):
        wf = WorkflowBuilder("cwl_test_array_step_input")
        wf.input("inp1", Optional[str])
        wf.input("inp2", Optional[str])

        wf.step(
            "print",
            ArrayTestTool(
                inps=[
                    If(IsDefined(wf.inp1), wf.inp1, "default1"),
                    If(IsDefined(wf.inp2), wf.inp2 + "_suffix", ""),
                ]
            ),
        ),

        wf.output("out", source=wf.print)

        ret, _, _ = wf.translate("wdl", to_console=False, allow_empty_container=True)

        expected = """\
version development

import "tools/ArrayStepTool.wdl" as A

workflow cwl_test_array_step_input {
  input {
    String? inp1
    String? inp2
  }
  call A.ArrayStepTool as print {
    input:
      inps=[if (defined(inp1)) then inp1 else "default1", if (defined(inp2)) then (inp2 + "_suffix") else ""]
  }
  output {
    Array[File] out = print.outs
  }
}"""

        self.assertEqual(expected, ret)

    def test_expression_for_output(self):
        ti = ToolInput("bam", Array(File))

        res = wdl.WdlTranslator.unwrap_expression(
            If(
                InputSelector("bam").length().equals(1),
                InputSelector("bam")[0],
                "generated",
            ),
            for_output=True,
            inputsdict={"bam": ti},
            string_environment=False,
        )

        self.assertEqual('if ((length(bam) == 1)) then bam[0] else "generated"', res)


class TestWdlWorkflowInputToOutputConnection(unittest.TestCase):
    def test_simple(self):
        w = WorkflowBuilder("wf")
        w.input("inp", str)
        w.output("out", source=w.inp)
        out, _, _ = w.translate("wdl", to_console=False)
        expected = """\
version development



workflow wf {
  input {
    String inp
  }
  output {
    String out = inp
  }
}"""
        self.assertEqual(expected, out)

    def test_with_int_default(self):
        w = WorkflowBuilder("wf")
        w.input("inp", int, default=0)
        w.output("out", source=w.inp)
        out, _, _ = w.translate("wdl", to_console=False)
        expected = """\
version development



workflow wf {
  input {
    Int? inp = 0
  }
  output {
    Int out = select_first([inp, 0])
  }
}"""
        self.assertEqual(expected, out)

    def test_with_str_default(self):
        w = WorkflowBuilder("wf")
        w.input("inp", str, default="hello")
        w.output("out", source=w.inp)
        out, _, _ = w.translate("wdl", to_console=False)
        expected = """\
version development



workflow wf {
  input {
    String? inp = "hello"
  }
  output {
    String out = select_first([inp, "hello"])
  }
}"""
        self.assertEqual(expected, out)


class TestWdlResourceOperators(unittest.TestCase):
    def test_1(self):
        tool_wdl = WdlTranslator.translate_tool_internal(
            OperatorResourcesTestTool(), with_resource_overrides=True
        ).get_string()
        lines = tool_wdl.splitlines(keepends=False)
        cpus = lines[-12].strip()
        time = lines[-9].strip()
        memory = lines[-8].strip()

        self.assertEqual("cpu: select_first([runtime_cpu, (2 * outputFiles), 1])", cpus)
        self.assertEqual(
            'memory: "~{select_first([runtime_memory, if ((size(inputFile, "MB") > 1024)) then 4 else 2, 4])}G"',
            memory,
        )
        self.assertEqual("duration: select_first([runtime_seconds, 60, 86400])", time)

    def test_base(self):
        tool_wdl = WdlTranslator.translate_tool_internal(
            EchoTestTool(), with_resource_overrides=True
        ).get_string()
        lines = tool_wdl.splitlines(keepends=False)
        # print(tool_wdl)
        cpus = lines[-12].strip()
        disks = lines[-11].strip()
        time = lines[-9].strip()
        memory = lines[-8].strip()

        self.assertEqual("cpu: select_first([runtime_cpu, 1])", cpus)

        self.assertEqual('memory: "~{select_first([runtime_memory, 4])}G"', memory)

        self.assertEqual("duration: select_first([runtime_seconds, 86400])", time)

        self.assertEqual(
            'disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"', disks
        )


class TestReadContentsOperator(unittest.TestCase):
    def test_read_contents_string(self):

        t = CommandToolBuilder(
            tool="test_readcontents",
            base_command=["echo", "1"],
            inputs=[ToolInput("inp", File)],
            outputs=[
                ToolOutput("out", String, glob=ReadContents(Stdout())),
                ToolOutput(
                    "out_json", String, selector=InputSelector("inp").read_json()["out"]
                ),
            ],
            container=None,
            version="-1",
        )

        t.translate("wdl", allow_empty_container=True)

    def test_read_contents_as_int(self):

        t = CommandToolBuilder(
            tool="test_readcontents",
            base_command=["echo", "1"],
            inputs=[],
            outputs=[ToolOutput("out", Float, glob=ReadContents(Stdout()).as_float())],
            container=None,
            version="-1",
        )

        t.translate("wdl", allow_empty_container=True)


class TestWDLNotNullOperator(unittest.TestCase):
    def test_workflow_string_not_null(self):
        w = WorkflowBuilder("wf")
        w.input("inp", Optional[str])
        w.output("out", source=w.inp.assert_not_null())

        wdltool = (
            w.translate("wdl", allow_empty_container=True, to_console=False)[0]
            .splitlines()[-3]
            .strip()
        )
        self.assertEqual("String out = select_first([inp])", wdltool)

    def test_commandtool_string(self):

        t = CommandToolBuilder(
            tool="id",
            base_command=None,
            inputs=[ToolInput("inp", Optional[str])],
            outputs=[
                ToolOutput("out", str, glob=InputSelector("inp").assert_not_null())
            ],
            version=None,
            container=None,
        )

        wdltool = (
            t.translate("wdl", allow_empty_container=True, to_console=False)
            .splitlines()[-3]
            .strip()
        )
        self.assertEqual("String out = select_first([inp])", wdltool)


class TestWdlWildcardSelector(unittest.TestCase):
    def test_regular_wildcard_selector(self):
        out = ToolOutput("out", Array(File), selector=WildcardSelector("*.txt"))
        translated_out = WdlTranslator.translate_tool_outputs([out], {}, out)[0]
        self.assertEqual('Array[File] out = glob("*.txt")', translated_out.get_string())

    def test_regular_wildcard_selector_single(self):
        out = ToolOutput(
            "out", File, selector=WildcardSelector("*.txt", select_first=True)
        )
        translated_out = WdlTranslator.translate_tool_outputs([out], {}, out)[0]
        self.assertEqual('File out = glob("*.txt")[0]', translated_out.get_string())

    def test_regular_wildcard_selector_single_warning(self):
        out = ToolOutput("out", File, selector=WildcardSelector("*.txt"))
        translated_out = WdlTranslator.translate_tool_outputs([out], {}, out)[0]
        self.assertEqual('File out = glob("*.txt")[0]', translated_out.get_string())

    def test_regular_wildcard_selector_single_optional(self):
        out = ToolOutput(
            "out",
            File(optional=True),
            selector=WildcardSelector("*.txt", select_first=True),
        )
        translated_out = WdlTranslator.translate_tool_outputs([out], {}, out)[0]
        self.assertEqual(
            'File? out = if length(glob("*.txt")) > 0 then glob("*.txt")[0] else None',
            translated_out.get_string(),
        )


class TestUnionType(unittest.TestCase):
    def test_lots_of_files(self):
        class TextFile(File):
            def name(self):
                return "TextFile"

        uniontype = UnionType(File, TextFile)
        self.assertEqual("File", uniontype.wdl().get_string())

    def test_file_int_fail(self):
        uniontype = UnionType(File, int)
        self.assertRaises(Exception, uniontype.wdl)


class TestForEachSelectors(unittest.TestCase):
    def test_minimal(self):
        with TemporaryDirectory() as tmpdir:
            TestForEach().translate("wdl", to_disk=True, export_path=tmpdir)
        w, _ = WdlTranslator.translate_workflow(TestForEach())
        expected = """\
version development

import "tools/EchoTestTool_TEST.wdl" as E

workflow TestForEach {
  input {
    Array[String] inp
  }
  scatter (idx in inp) {
     call E.EchoTestTool as print {
      input:
        inp=(idx + "-hello")
    }
  }
  output {
    Array[File] out = print.out
  }
}"""
        self.assertEqual(expected.strip(), w.get_string().strip())


t = CommandToolBuilder(
    tool="test_readcontents",
    base_command=["echo", "1"],
    inputs=[ToolInput("inp", File)],
    outputs=[
        ToolOutput(
            "out_json", String, selector=InputSelector("inp").read_json()["out"]
        ),
    ],
    container=None,
    version="-1",
)
