import unittest
from typing import List, Dict, Any

import wdlgen
from janis_core.utils.scatter import ScatterDescription, ScatterMethod, ScatterMethods

import janis_core.translations.wdl as wdl
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
    Filename,
    WildcardSelector,
    ToolArgument,
    Boolean,
    Int,
)
from janis_core.tests.testtools import SingleTestTool
from janis_core.translations import WdlTranslator
from janis_core.operators.selectors import CpuSelector, StringFormatter


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


class TestTool(CommandTool):
    @staticmethod
    def tool():
        return "TestTranslationtool"

    @staticmethod
    def base_command():
        return "echo"

    def inputs(self) -> List[ToolInput]:
        return [ToolInput("testtool", String()), ToolInput("arrayInp", Array(File()))]

    def arguments(self) -> List[ToolArgument]:
        return [ToolArgument(StringFormatter('test:\\t:escaped:\\n:characters"'))]

    def outputs(self) -> List[ToolOutput]:
        return [ToolOutput("std", Stdout())]

    def cpus(self, hints: Dict[str, Any]):
        return 2

    def memory(self, hints: Dict[str, Any]):
        return 2

    def friendly_name(self) -> str:
        return "Tool for testing translation"

    @staticmethod
    def container():
        return "ubuntu:latest"

    @staticmethod
    def version():
        return None

    def env_vars(self):
        return {"test1": InputSelector("testtool")}


class TestToolWithSecondaryOutput(TestTool):
    def outputs(self):
        return [
            ToolOutput(
                "out",
                TestTypeWithNonEscapedSecondary(),
                glob=InputSelector("testtool") + "/out",
            )
        ]


class TestTypeWithSecondary(File):
    @staticmethod
    def secondary_files():
        return ["^.txt"]


class TestTypeWithNonEscapedSecondary(File):
    @staticmethod
    def secondary_files():
        return [".txt"]


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

        outp = wdl.translate_step_node(stp2, stp2.id(), {}, set())
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
        import uuid

        fn = Filename(guid=str(uuid.uuid4()))
        self.assertEqual(
            fn.generated_filename(),
            wdl.WdlTranslator.unwrap_expression(fn, None, string_environment=True),
        )

    def test_input_value_filename_nostringenv(self):
        import uuid

        fn = Filename(guid=str(uuid.uuid4()))
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

    def test_input_value_cpuselect_stringenv(self):
        # CpuSelector relies on their being a runtime_cpu attribute,
        # this test will assume it's present, and '' will test ensure
        # that it is actually present
        ti = {"runtime_cpu": ToolInput("runtime_cpu", Int(), default=1)}
        inp = CpuSelector()
        self.assertEqual(
            "~{select_first([runtime_cpu, 1])}",
            wdl.WdlTranslator.unwrap_expression(inp, ti, string_environment=True),
        )

    def test_input_value_cpuselect_nostringenv(self):
        # CpuSelector relies on their being a runtime_cpu attribute,
        # this test will assume it's present, and '' will test ensure
        # that it is actually present

        ti = {"runtime_cpu": ToolInput("runtime_cpu", Int(), default=1)}
        inp = CpuSelector()
        self.assertEqual(
            "select_first([runtime_cpu, 1])",
            wdl.WdlTranslator.unwrap_expression(inp, ti, string_environment=False),
        )

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

    def test_string_formatter_optional_inpselect_no_default(self):
        # will throw
        ti = {"ti": ToolInput("ti", String(optional=True))}
        b = StringFormatter("{place} michael", place=InputSelector("ti"))
        self.assertRaises(Exception, wdl.WdlTranslator.unwrap_expression, b, ti)

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
        os = WdlTranslator.translate_tool_outputs([toolout], inmap, toolid=tool.id())

        self.assertEqual("out", os[0].name)
        self.assertEqual('(testtool + "/out")', os[0].expression)

        self.assertEqual("out_txt", os[1].name)
        self.assertEqual('(testtool + "/out") + ".txt"', os[1].expression)


class TestWdlToolInputGeneration(unittest.TestCase):
    def test_nodefault_nooptional_position(self):
        ti = ToolInput("tag", String(), position=0)
        resp = wdl.translate_command_input(ti)
        self.assertEqual("~{tag}", resp.get_string())

    def test_nodefault_nooptional_prefix_sep(self):
        ti = ToolInput("tag", String(), prefix="--amazing")
        resp = wdl.translate_command_input(ti)
        self.assertEqual("--amazing ~{tag}", resp.get_string())

    def test_nodefault_nooptional_prefix_nosep(self):
        ti = ToolInput(
            "tag", String(), prefix="--amazing=", separate_value_from_prefix=False
        )
        resp = wdl.translate_command_input(ti)
        self.assertEqual("--amazing=~{tag}", resp.get_string())

    def test_nodefault_optional_position(self):
        ti = ToolInput("tag", String(optional=True), position=0)
        resp = wdl.translate_command_input(ti)
        self.assertEqual("~{tag}", resp.get_string())

    def test_nodefault_optional_prefix_sep(self):
        ti = ToolInput("tag", String(optional=True), prefix="--amazing")
        resp = wdl.translate_command_input(ti)
        self.assertEqual(
            '~{if defined(tag) then ("--amazing " +  \'"\' + tag + \'"\') else ""}',
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
            '~{if defined(tag) then (\'"\' + "--amazing=" + tag + \'"\') else ""}',
            resp.get_string(),
        )

    def test_default_nooptional_position(self):
        # this will get turned into an optional
        ti = ToolInput("tag", String(), position=0, default="defval")
        resp = wdl.translate_command_input(ti)
        self.assertEqual('~{select_first([tag, "defval"])}', resp.get_string())

    def test_default_nooptional_prefix_sep(self):
        ti = ToolInput("tag", String(), prefix="--amazing", default="defval")
        resp = wdl.translate_command_input(ti)
        self.assertEqual(
            '--amazing ~{select_first([tag, "defval"])}', resp.get_string()
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
            '--amazing=~{select_first([tag, "defval"])}', resp.get_string()
        )

    def test_default_optional_position(self):
        ti = ToolInput("tag", String(optional=True), position=0, default="defval")
        resp = wdl.translate_command_input(ti)
        self.assertEqual('~{select_first([tag, "defval"])}', resp.get_string())

    def test_default_optional_prefix_sep(self):
        ti = ToolInput(
            "tag", String(optional=True), prefix="--amazing", default="defval"
        )
        resp = wdl.translate_command_input(ti)
        self.assertEqual(
            '~{if defined(select_first([tag, "defval"])) then ("--amazing " +  \'"\' + select_first([tag, "defval"]) + \'"\') else ""}',
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
            '~{if defined(select_first([tag, "defval"])) then (\'"\' + "--amazing=" + select_first([tag, "defval"]) + \'"\') else ""}',
            resp.get_string(),
        )

    def test_bind_boolean_as_default(self):
        ti = ToolInput("tag", Boolean(optional=True), prefix="--amazing", default=True)
        resp = wdl.translate_command_input(ti).get_string()
        self.assertEqual('~{true="--amazing" false="" select_first([tag, true])}', resp)


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
        resources = WdlTranslator.build_resources_input(tool.wrapped_in_wf(), {})
        self.assertEqual(
            2, resources["TestTranslationtoolWf.testtranslationtool_runtime_cpu"]
        )

    def test_max_cores(self):
        tool = TestTool()
        resources = WdlTranslator.build_resources_input(
            tool.wrapped_in_wf(), {}, max_cores=1
        )
        self.assertEqual(
            1, resources["TestTranslationtoolWf.testtranslationtool_runtime_cpu"]
        )

    def test_memory(self):
        tool = TestTool()
        resources = WdlTranslator.build_resources_input(tool.wrapped_in_wf(), {})
        self.assertEqual(
            2, resources["TestTranslationtoolWf.testtranslationtool_runtime_memory"]
        )

    def test_max_memory(self):
        tool = TestTool()
        resources = WdlTranslator.build_resources_input(
            tool.wrapped_in_wf(), {}, max_mem=1
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
            SingleTestTool(inputs=w.inp, input2=w.inp2),
            scatter=ScatterDescription(fields=["inputs"], method=ScatterMethods.dot),
        )

        outp = wdl.translate_step_node(step, "A.SingleTestTool", {}, {"inp", "inp2"})
        expected = """\
scatter (i in inp) {
   call A.SingleTestTool as dotTool {
    input:
      inputs=i,
      input2=inp2
  }
}"""
        self.assertEqual(expected, outp.get_string(indent=0))

    def test_scatter_single_no_description(self):
        w = WorkflowBuilder("sbmf")
        w.input("inp", Array(str))
        w.input("inp2", str)

        step = w.step(
            "dotTool", SingleTestTool(inputs=w.inp, input2=w.inp2), scatter="inputs"
        )

        outp = wdl.translate_step_node(step, "A.SingleTestTool", {}, {"inp", "inp2"})
        expected = """\
scatter (i in inp) {
   call A.SingleTestTool as dotTool {
    input:
      inputs=i,
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
            SingleTestTool(inputs=w.inp, input2=w.inp2),
            scatter=ScatterDescription(
                fields=["inputs", "input2"], method=ScatterMethods.dot
            ),
        )

        outp = wdl.translate_step_node(step, "A.SingleTestTool", {}, {"inp", "inp2"})
        expected = """\
scatter (Q in zip(inp, inp2)) {
   call A.SingleTestTool as dotTool {
    input:
      inputs=Q.left,
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
            SingleTestTool(inputs=w.inp, input2=w.inp2, input3=w.inp3),
            scatter=ScatterDescription(
                fields=["inputs", "input2", "input3"], method=ScatterMethods.dot
            ),
        )

        outp = wdl.translate_step_node(
            step, "A.SingleTestTool", {}, {"inp", "inp2", "inp3"}
        )
        expected = """\
scatter (Q in zip(inp, zip(inp2, inp3))) {
   call A.SingleTestTool as dotTool {
    input:
      inputs=Q.left,
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
            SingleTestTool(inputs=w.inp, input2=w.inp2, input3=w.inp3, input4=w.inp4),
            scatter=ScatterDescription(
                fields=["inputs", "input2", "input3", "input4"],
                method=ScatterMethods.dot,
            ),
        )

        outp = wdl.translate_step_node(
            step, "A.SingleTestTool", {}, {"inp", "inp2", "inp3", "inp4"}
        )
        expected = """\
scatter (Q in zip(inp, zip(inp2, zip(inp3, inp4)))) {
   call A.SingleTestTool as dotTool {
    input:
      inputs=Q.left,
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
                fields=["input1", "input2"], method=ScatterMethods.dot
            ),
        )

        outp = wdl.translate_step_node(step, "A.SingleTestTool", {}, {"inp", "inp2"})
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
        w.step("echo", SingleTestTool(inputs=w.inp))
        w.step("echo_2", SingleTestTool(inputs=w.inp))

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
    String echo_runtime_disks
    Int? echo_2_runtime_memory
    Int? echo_2_runtime_cpu
    String echo_2_runtime_disks
  }
  call T.TestStepTool as echo {
    input:
      inputs=inp,
      runtime_memory=echo_runtime_memory,
      runtime_cpu=echo_runtime_cpu,
      runtime_disks=echo_runtime_disks
  }
  call T.TestStepTool as echo_2 {
    input:
      inputs=inp,
      runtime_memory=echo_2_runtime_memory,
      runtime_cpu=echo_2_runtime_cpu,
      runtime_disks=echo_2_runtime_disks
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
