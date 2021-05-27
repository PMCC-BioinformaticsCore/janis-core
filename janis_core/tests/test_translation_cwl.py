import unittest
from typing import List, Dict, Any, Optional

from janis_core.tool.commandtool import ToolArgument

from janis_core.operators.logical import If, IsDefined
from janis_core.operators.standard import ReadContents, FilterNullOperator

from janis_core.tests.testtools import (
    SingleTestTool,
    ArrayTestTool,
    TestTool,
    TestToolV2,
    TestTypeWithSecondary,
    TestWorkflowWithStepInputExpression,
    EchoTestTool,
    FilenameGeneratedTool,
    OperatorResourcesTestTool,
    TestForEach,
)

from janis_core.deps import cwlgen

import janis_core.translations.cwl as cwl
from janis_core import (
    WorkflowBuilder,
    ToolInput,
    String,
    InputSelector,
    Array,
    WildcardSelector,
    StringFormatter,
    CommandToolBuilder,
    ToolOutput,
    DataType,
    Float,
)
from janis_core.tool.documentation import InputDocumentation
from janis_core.translations import CwlTranslator
from janis_core.types import CpuSelector, MemorySelector, Stdout, UnionType, File
from janis_core.workflow.workflow import InputNode


class TestTypeWithAlternateAndSecondary(File):
    def __init__(self, optional=False):
        super().__init__(
            optional=optional, extension=".txt", alternate_extensions={".text"}
        )

    @staticmethod
    def secondary_files():
        return ["^.file"]


class TestCwlTypesConversion(unittest.TestCase):
    pass


class TestCwlMisc(unittest.TestCase):
    def test_str_tool(self):
        t = TestTool()
        actual = t.translate("cwl", to_console=False)
        self.maxDiff = None
        self.assertEqual(cwl_testtool, actual)


class TestCwlTranslatorOverrides(unittest.TestCase):
    def setUp(self):
        self.translator = CwlTranslator()

    def test_stringify_WorkflowBuilder(self):
        cwlobj = cwlgen.Workflow(
            id="wid", cwlVersion="v1.2", inputs={}, outputs={}, steps={}
        )
        expected = """\
#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.2

inputs: {}

outputs: {}

steps: {}
id: wid
"""
        self.assertEqual(
            expected, self.translator.stringify_translated_workflow(cwlobj)
        )

    def test_stringify_tool(self):
        cwlobj = cwlgen.CommandLineTool(
            id="tid", inputs={}, outputs={}, cwlVersion="v1.2"
        )
        expected = """\
#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.2

inputs: {}

outputs: {}
id: tid
"""

        self.assertEqual(expected, self.translator.stringify_translated_tool(cwlobj))

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
        cwltoolinput = cwl.translate_tool_input(t, None, None).save()
        self.assertDictEqual(
            {
                "id": "filesA",
                "label": "filesA",
                "type": {"items": "string", "type": "array"},
                "inputBinding": {"prefix": "-A", "position": 1},
            },
            cwltoolinput,
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
        cwltoolinput = cwl.translate_tool_input(t, None, None)
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
            cwltoolinput.save(),
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
        cwltoolinput = cwl.translate_tool_input(t, None, None)
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
            cwltoolinput.save(),
        )

    def test_optional_array_prefixes(self):
        t = ToolInput(
            "filesD",
            Array(String(), optional=True),
            prefix="-D",
            prefix_applies_to_all_elements=True,
        )
        cwltoolinput = cwl.translate_tool_input(t, None, None)

        self.assertDictEqual(
            {
                "id": "filesD",
                "label": "filesD",
                "inputBinding": {},
                "type": [
                    {
                        "inputBinding": {"prefix": "-D"},
                        "items": "string",
                        "type": "array",
                    },
                    "null",
                ],
            },
            dict(cwltoolinput.save()),
        )


class TestCwlSelectorsAndGenerators(unittest.TestCase):
    def test_input_selector_base(self):
        input_sel = InputSelector("random")
        self.assertEqual(
            "$(inputs.random)",
            cwl.translate_input_selector(
                input_sel,
                code_environment=False,
                inputs_dict={},
                skip_inputs_lookup=True,
            ),
        )

    def test_input_selector_base_codeenv(self):
        input_sel = InputSelector("random")
        self.assertEqual(
            "inputs.random",
            cwl.translate_input_selector(
                input_sel,
                code_environment=True,
                inputs_dict={},
                skip_inputs_lookup=True,
            ),
        )

    def test_input_value_none_codeenv(self):
        self.assertEqual(
            "null", cwl.CwlTranslator.unwrap_expression(None, code_environment=True)
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
            "42", cwl.CwlTranslator.unwrap_expression(42, code_environment=True)
        )

    def test_input_value_int_nocodeenv(self):
        self.assertEqual(
            "42", cwl.CwlTranslator.unwrap_expression(42, code_environment=False)
        )

    def test_alias_selector(self):
        w = WorkflowBuilder("wf")
        w.input("inp", str)
        w.step("echo", EchoTestTool(inp=w.inp.as_type(str)))
        w.output("out", source=w.echo.out)
        sn: List[cwlgen.WorkflowStep] = cwl.translate_step_node(
            w.step_nodes["echo"], inputs_dict={"inp": ToolInput("inp", str)}
        )
        self.assertEqual("inp", sn[0].in_[0].source)

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
            cwl.CwlTranslator.unwrap_expression(
                inp,
                code_environment=True,
                inputs_dict={"threads": ToolInput("threads", int)},
            ),
        )

    def test_input_value_inpselect_nocodeenv(self):
        inp = InputSelector("threads")
        self.assertEqual(
            "$(inputs.threads)",
            cwl.CwlTranslator.unwrap_expression(
                inp,
                code_environment=False,
                inputs_dict={"threads": ToolInput("threads", int)},
            ),
        )

    def test_input_value_removing_extension(self):
        clt = CommandToolBuilder(
            tool="dev",
            base_command="echo",
            inputs=[ToolInput("inp", File(extension=".txt"))],
            outputs=[ToolOutput("out", str, selector=ReadContents(Stdout()))],
            version="v1.0",
            container="ubuntu",
        )

        arg = cwl.translate_tool_argument(
            ToolArgument(InputSelector("inp", remove_file_extension=True) + ".bam"),
            clt,
            inputs_dict={t.id(): t for t in clt.inputs()},
        )

        self.assertEqual(
            '$((inputs.inp.basename.replace(/.txt$/, "") + ".bam"))', arg.valueFrom
        )

    # def test_input_value_wildcard(self):
    #     self.assertRaises(
    #         Exception, cwl.CwlTranslator.unwrap_expression, value=WildcardSelector("*")
    #     )

    # def test_input_value_cpuselect_codeenv(self):
    #     inp = CpuSelector()
    #     self.assertEqual(
    #         "inputs.runtime_cpu",
    #         cwl.CwlTranslator.unwrap_expression(inp, code_environment=True),
    #     )
    #
    # def test_input_value_cpuselect_nocodeenv(self):
    #     inp = CpuSelector()
    #     self.assertEqual(
    #         "$(inputs.runtime_cpu)",
    #         cwl.CwlTranslator.unwrap_expression(inp, code_environment=False),
    #     )

    # def test_input_value_memselect_codeenv(self):
    #     inp = MemorySelector()
    #     self.assertEqual(
    #         "inputs.runtime_memory",
    #         cwl.CwlTranslator.unwrap_expression(inp, code_environment=True),
    #     )
    #
    # def test_input_value_memselect_nocodeenv(self):
    #     inp = MemorySelector()
    #     self.assertEqual(
    #         "$(inputs.runtime_memory)",
    #         cwl.CwlTranslator.unwrap_expression(inp, code_environment=False),
    #     )

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
        res = cwl.CwlTranslator.unwrap_expression(b, code_environment=False)
        self.assertEqual('$("there\'s {one} arg".replace(/\{one\}/g, "a string"))', res)

    def test_string_formatter_one_input_selector_param(self):
        b = StringFormatter("an input {arg}", arg=InputSelector("random_input"))
        res = cwl.CwlTranslator.unwrap_expression(
            b,
            code_environment=False,
            inputs_dict={"random_input": ToolInput("random_input", str)},
        )
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
        inputs_dict = {
            "tumorInputName": ToolInput("tumorInputName", str),
            "normalInputName": ToolInput("normalInputName", str),
        }
        res = cwl.CwlTranslator.unwrap_expression(
            b, code_environment=False, inputs_dict=inputs_dict
        )
        self.assertEqual(
            '$("{tumorName}:{normalName}".replace(/\{tumorName\}/g, inputs.tumorInputName).replace(/\{normalName\}/g, inputs.normalInputName))',
            res,
        )

    def test_escaped_characters(self):
        trans = cwl.CwlTranslator
        translated = trans.translate_tool_internal(TestTool())
        arg: cwlgen.CommandLineBinding = translated.arguments[0]
        self.assertEqual('test:\\\\t:escaped:\\\\n:characters\\"', arg.valueFrom)


class TestCwlEnvVar(unittest.TestCase):
    def test_environment1(self):
        t = CwlTranslator().translate_tool_internal(tool=TestTool())
        envvar: cwlgen.EnvVarRequirement = [
            t for t in t.requirements if t.class_ == "EnvVarRequirement"
        ][0]
        envdef: cwlgen.EnvironmentDef = envvar.envDef[0]
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
        tinp = cwl.translate_workflow_input(inp, None)

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
        tinp = cwl.translate_workflow_input(inp, None)

        self.assertEqual("File", tinp.type)
        self.assertEqual("^.txt", tinp.secondaryFiles[0].pattern)

    def test_array_secondary_file_translation(self):
        inp = InputNode(
            None,
            identifier="testIdentifier",
            datatype=Array(TestTypeWithSecondary()),
            default=None,
            value=None,
        )
        tinp = cwl.translate_workflow_input(inp, None)
        self.assertIsInstance(tinp.type, cwlgen.CommandInputArraySchema)
        self.assertEqual("File", tinp.type.items)
        self.assertEqual("^.txt", tinp.secondaryFiles[0].pattern)


class TestCwlOutputGeneration(unittest.TestCase):
    def test_stdout_no_outputbinding(self):
        out = cwl.translate_tool_output(ToolOutput("out", Stdout), {}, tool=None).save()
        self.assertDictEqual({"id": "out", "label": "out", "type": "stdout"}, out)

    def test_localised_out(self):

        inps = {"inp": ToolInput("inp", File, position=1, localise_file=True)}
        out = ToolOutput("out", File, selector=InputSelector("inp"))

        cwlout = cwl.translate_tool_output(
            out, inps, environment="dev-test_localised_out", tool=None
        )
        ob: cwlgen.CommandOutputBinding = cwlout.outputBinding
        self.assertEqual("$(inputs.inp.basename)", ob.glob)


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
        self.assertDictEqual({"inpId": None}, self.translator.build_inputs_file(wf))

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
        w.step("stp1", ArrayTestTool(inps=w.inp1))

        c, _, _ = CwlTranslator().translate(
            w, to_console=False, allow_empty_container=True
        )
        self.assertEqual(cwl_multiinput, c)


class TestPackedWorkflow(unittest.TestCase):
    def test_simple(self):
        w = WorkflowBuilder("test_add_single_to_array_edge")
        w.step("ech", SingleTestTool(input1="Hello"), doc="Print 'Hello'")
        c = CwlTranslator.translate_workflow_to_all_in_one(
            w, allow_empty_container=True
        )
        print(CwlTranslator.stringify_translated_workflow(c))


class TestContainerOverride(unittest.TestCase):
    def test_tool_dict_override(self):
        import ruamel.yaml

        expected_container = "container/override"

        tool = SingleTestTool()
        d = ruamel.yaml.load(
            tool.translate(
                "cwl",
                to_console=False,
                container_override={tool.id(): expected_container},
            ),
            Loader=ruamel.yaml.Loader,
        )

        received_container = [
            req.get("dockerPull")
            for req in d.get("requirements")
            if req["class"] == "DockerRequirement"
        ][0]

        self.assertEqual(expected_container, received_container)

    def test_tool_string_override(self):
        import ruamel.yaml

        expected_container = "container/override"

        tool = SingleTestTool()
        d = ruamel.yaml.load(
            tool.translate(
                "cwl", to_console=False, container_override=expected_container
            ),
            Loader=ruamel.yaml.Loader,
        )
        received_container = [
            req.get("dockerPull")
            for req in d.get("requirements")
            if req["class"] == "DockerRequirement"
        ][0]

        self.assertEqual(expected_container, received_container)


class TestCWLCompleteOperators(unittest.TestCase):
    def test_step_input(self):
        self.maxDiff = None

        ret, _, _ = TestWorkflowWithStepInputExpression().translate(
            "cwl", to_console=False
        )
        self.assertEqual(cwl_stepinput, ret)

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

        ret, _, _ = wf.translate("cwl", allow_empty_container=True, to_console=False)
        self.maxDiff = None
        self.assertEqual(cwl_arraystepinput, ret)


class TestCreateFilesAndDirectories(unittest.TestCase):

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
        req = CwlTranslator.build_initial_workdir_from_tool(command).listing

        self.assertEqual(1, len(req))
        self.assertEqual(
            '$({ class: "Directory", basename: "test-directory", listing: [] })', req[0]
        )

    def test_create_single_directory_from_selector(self):
        command = CommandToolBuilder(
            **self.initial_params, directories_to_create=InputSelector("name")
        )
        req = CwlTranslator.build_initial_workdir_from_tool(command).listing
        self.assertEqual(1, len(req))
        self.assertEqual(
            '$({ class: "Directory", basename: inputs.name, listing: [] })', req[0]
        )

    def test_create_single_directory_from_operator(self):
        command = CommandToolBuilder(
            **self.initial_params, directories_to_create=InputSelector("name") + "-out"
        )
        req = CwlTranslator.build_initial_workdir_from_tool(command).listing
        self.assertEqual(1, len(req))
        self.assertEqual(
            '$({ class: "Directory", basename: (inputs.name + "-out"), listing: [] })',
            req[0],
        )

    def test_create_single_file_from_operator(self):
        command = CommandToolBuilder(
            **self.initial_params,
            files_to_create=[("my-path.txt", InputSelector("inp").contents())],
        )
        req = CwlTranslator.build_initial_workdir_from_tool(command).listing
        self.assertEqual(1, len(req))
        self.assertEqual("my-path.txt", req[0].entryname)
        self.assertEqual("$(inputs.inp.contents)", req[0].entry)

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
        req = CwlTranslator.build_initial_workdir_from_tool(command).listing
        self.assertEqual(1, len(req))
        self.assertIsInstance(req[0], cwlgen.Dirent)
        self.assertEqual(
            '$("{name}.txt".replace(/\{name\}/g, inputs.name))', req[0].entryname
        )
        self.assertEqual("this is contents", req[0].entry)


class WorkflowCwlInputDefaultOperator(unittest.TestCase):
    def test_string_formatter(self):
        wf = WorkflowBuilder("wf")
        wf.input("sampleName", str)
        wf.input("platform", str)

        wf.input(
            "readGroupHeaderLine",
            String(optional=True),
            default=StringFormatter(
                "@RG\\tID:{name}\\tSM:{name}\\tLB:{name}\\tPL:{pl}",
                name=InputSelector("sampleName"),
                pl=InputSelector("platform"),
            ),
        )
        wf.step("print", EchoTestTool(inp=wf.readGroupHeaderLine))
        wf.output("out", source=wf.print)
        d, _ = cwl.CwlTranslator.translate_workflow(
            wf, with_container=False, allow_empty_container=True
        )
        stepinputs = d.save()["steps"][0]["in"]
        self.assertEqual(4, len(stepinputs))
        expression = stepinputs[-1]["valueFrom"]
        expected = (
            "$((inputs._print_inp_readGroupHeaderLine != null) "
            "? inputs._print_inp_readGroupHeaderLine "
            ': "@RG\\\\tID:{name}\\\\tSM:{name}\\\\tLB:{name}\\\\tPL:{pl}".replace(/\\{name\\}/g, inputs._print_inp_sampleName).replace(/\\{pl\\}/g, inputs._print_inp_platform))'
        )
        self.assertEqual(expected, expression)

    def test_string_formatter_stepinput(self):
        wf = WorkflowBuilder("wf")
        wf.input("sampleName", str)
        wf.input("platform", str)

        wf.step(
            "print",
            EchoTestTool(
                inp=StringFormatter(
                    "@RG\\tID:{name}\\tSM:{name}\\tLB:{name}\\tPL:{pl}",
                    name=wf.sampleName,
                    pl=wf.platform,
                )
            ),
        )
        wf.output("out", source=wf.print)
        d, _ = cwl.CwlTranslator.translate_workflow(
            wf, with_container=False, allow_empty_container=True
        )
        stepinputs = d.save()["steps"][0]["in"]
        self.assertEqual(3, len(stepinputs))
        expression = stepinputs[-1]["valueFrom"]
        expected = '$("@RG\\\\tID:{name}\\\\tSM:{name}\\\\tLB:{name}\\\\tPL:{pl}".replace(/\\{name\\}/g, inputs._print_inp_sampleName).replace(/\\{pl\\}/g, inputs._print_inp_platform))'
        self.assertEqual(expected, expression)


class TestCWLFilenameGeneration(unittest.TestCase):
    def test_1(self):
        tool = FilenameGeneratedTool()
        inputsdict = {t.id(): t for t in tool.inputs()}
        mapped = [cwl.translate_tool_input(i, inputsdict, tool) for i in tool.inputs()]
        expressions = [
            mapped[i].save()["inputBinding"]["valueFrom"] for i in range(4, len(mapped))
        ]
        self.assertEqual("$(inputs.inp)", expressions[0])
        self.assertEqual(
            '$(inputs.inpOptional ? inputs.inpOptional : "generated")', expressions[1]
        )
        self.assertEqual(
            '$(inputs.fileInp.basename.replace(/.txt$/, "")).transformed.fnp',
            expressions[2],
        )
        self.assertEqual(
            '$(inputs.fileInpOptional ? inputs.fileInpOptional.basename.replace(/.txt$/, "") : "generated").optional.txt',
            expressions[3],
        )


class TestCWLRunRefs(unittest.TestCase):
    def test_two_similar_tools(self):
        w = WorkflowBuilder("testTwoToolsWithSameId")

        w.input("inp", str)
        w.step("stp1", TestTool(testtool=w.inp))
        w.step("stp2", TestToolV2(testtool=w.inp))

        wf_cwl, _ = CwlTranslator.translate_workflow(w)
        stps = {stp.id: stp for stp in wf_cwl.steps}

        self.assertEqual("tools/TestTranslationtool.cwl", stps["stp1"].run)
        self.assertEqual("tools/TestTranslationtool_v0_0_2.cwl", stps["stp2"].run)


class TestCwlResourceOperators(unittest.TestCase):
    def test_1(self):
        tool_cwl = CwlTranslator.translate_tool_internal(
            OperatorResourcesTestTool(), with_resource_overrides=True
        )
        resourcereq = [
            r for r in tool_cwl.requirements if r.class_ == "ResourceRequirement"
        ][0]
        self.assertEqual(
            "$([inputs.runtime_cpu, (2 * inputs.outputFiles), 1].filter(function (inner) { return inner != null })[0])",
            resourcereq.coresMin,
        )
        self.assertEqual(
            "$(Math.round((953.674 * [inputs.runtime_memory, ((inputs.inputFile.size / 1048576) > 1024) ? 4 : 2, 4].filter(function (inner) { return inner != null })[0])))",
            resourcereq.ramMin,
        )


class TestReadContentsOperator(unittest.TestCase):
    def test_read_contents_string(self):

        t = CommandToolBuilder(
            tool="test_readcontents",
            base_command=["echo", "1"],
            inputs=[],
            outputs=[ToolOutput("out", String, glob=ReadContents(Stdout()))],
            container=None,
            version="-1",
        )

        translated = CwlTranslator.translate_tool_internal(
            t, allow_empty_container=True
        )
        self.assertTrue(translated.outputs[0].outputBinding.loadContents)

    def test_read_contents_as_int(self):

        t = CommandToolBuilder(
            tool="test_readcontents",
            base_command=["echo", "1"],
            inputs=[],
            outputs=[ToolOutput("out", Float, glob=ReadContents(Stdout()).as_float())],
            container=None,
            version="-1",
        )
        translated = CwlTranslator.translate_tool_internal(
            t, allow_empty_container=True
        )
        self.assertTrue(translated.outputs[0].outputBinding.loadContents)
        self.assertEqual("float", translated.outputs[0].type)


class TestCWLNotNullOperator(unittest.TestCase):
    def test_workflow_string_not_null(self):
        w = WorkflowBuilder("wf")
        w.input("inp", Optional[str])
        w.output("out", source=w.inp.assert_not_null())

        cwltool = w.translate("cwl", allow_empty_container=True, to_console=False)[0]
        print(cwltool)

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

        cwltool = t.translate("cwl", allow_empty_container=True, to_console=False)
        print(cwltool)


class TestCwlScatterExpression(unittest.TestCase):
    def test_filter_null(self):
        T = CommandToolBuilder(
            tool="testsingleinput",
            base_command="echo",
            inputs=[ToolInput("inp", str, position=0)],
            outputs=[ToolOutput("out", Stdout)],
            version="v1",
            container=None,
        )
        w = WorkflowBuilder("wf")
        w.input("inp", Array(Optional[str], optional=True))
        w.step("stp", T(inp=FilterNullOperator(w.inp)), scatter="inp")
        w.output("out", source=w.stp.out)

        w_cwl = cwl.CwlTranslator().translate_workflow(w, with_container=False)[0]
        self.assertEqual(2, len(w_cwl.steps))
        self.assertEqual(
            "_evaluate_prescatter-stp-inp/out", w_cwl.steps[1].in_[0].source
        )


class TestWorkflowOutputExpression(unittest.TestCase):
    def test_read_contents(self):
        w = WorkflowBuilder("wf")
        w.input("inp", str)
        w.step("stp", EchoTestTool(inp=w.inp))
        w.output("out", source=w.stp.out.contents())

        w_cwl = cwl.CwlTranslator().translate_workflow(w, with_container=False)[0]

        self.assertEqual(2, len(w_cwl.steps))
        self.assertEqual(
            "${return {out: inputs._stpout.contents }}", w_cwl.steps[1].run.expression
        )
        self.assertTrue(w_cwl.steps[1].run.inputs[0].loadContents)


class TestCwlUnionType(unittest.TestCase):
    def test_file_file(self):
        utype = UnionType(File, File)
        cwl_utype = utype.cwl_type()
        self.assertEqual("File", cwl_utype)

    def test_file_int_str(self):
        utype = UnionType(int, File, File, str)
        cwl_utype = sorted(utype.cwl_type())
        self.assertListEqual(["File", "int", "string"], cwl_utype)


class TestCWLWhen(unittest.TestCase):
    def test_basic(self):
        w = WorkflowBuilder("my_conditional_workflow")

        w.input("inp", String(optional=True))

        w.step(
            "print_if_has_value",
            TestTool(testtool=w.inp),
            # only print if the input "inp" is defined.
            when=IsDefined(w.inp),
        )

        w.output("out", source=w.print_if_has_value)

        inputs_dict = {"inp": ToolInput("inp", str)}

        c = cwl.translate_step_node(w.print_if_has_value, inputs_dict=inputs_dict)[0]

        self.assertEqual("$((inputs.__when_inp != null))", c.when)
        extra_input: cwlgen.WorkflowStepInput = c.in_[-1]
        self.assertEqual("__when_inp", extra_input.id)


class TestForEachSelectors(unittest.TestCase):
    def test_minimal(self):
        tool = TestForEach()
        # tool.translate("cwl", export_path="~/Desktop/tmp", to_disk=True)
        w, _ = CwlTranslator.translate_workflow(tool)

        stp = w.steps[0]
        self.assertEqual("inp", stp.in_[0].source)
        self.assertEqual('$((inputs._idx + "-hello"))', stp.in_[1].valueFrom)


cwl_testtool = """\
#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.2
label: Tool for testing translation

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
- class: EnvVarRequirement
  envDef:
  - envName: test1
    envValue: $(inputs.testtool)
- class: DockerRequirement
  dockerPull: ubuntu:latest

inputs:
- id: testtool
  label: testtool
  type: string
- id: arrayInp
  label: arrayInp
  type:
  - type: array
    items: string
  - 'null'

outputs:
- id: std
  label: std
  type: stdout
stdout: _stdout
stderr: _stderr

baseCommand: echo
arguments:
- position: 0
  valueFrom: test:\\\\t:escaped:\\\\n:characters\\"

hints:
- class: ToolTimeLimit
  timelimit: |-
    $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
id: TestTranslationtool
"""


cwl_multiinput = """\
#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.2

requirements:
- class: InlineJavascriptRequirement
- class: StepInputExpressionRequirement
- class: MultipleInputFeatureRequirement

inputs:
- id: inp1
  type: string

outputs: []

steps:
- id: stp1
  in:
  - id: inps
    source:
    - inp1
    linkMerge: merge_nested
  run: tools/ArrayStepTool.cwl
  out:
  - id: outs
id: test_add_single_to_array_edge
"""

cwl_stepinput = """\
#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.2
label: 'TEST: WorkflowWithStepInputExpression'

requirements:
- class: InlineJavascriptRequirement
- class: StepInputExpressionRequirement

inputs:
- id: mystring
  type:
  - string
  - 'null'
- id: mystring_backup
  type:
  - string
  - 'null'

outputs:
- id: out
  type: File
  outputSource: print/out

steps:
- id: print
  in:
  - id: _print_inp_mystring
    source: mystring
  - id: _print_inp_mystringbackup
    source: mystring_backup
  - id: inp
    valueFrom: |-
      $((inputs._print_inp_mystring != null) ? inputs._print_inp_mystring : inputs._print_inp_mystringbackup)
  run: tools/EchoTestTool_TEST.cwl
  out:
  - id: out
id: TestWorkflowWithStepInputExpression
"""

cwl_arraystepinput = """\
#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.2

requirements:
- class: InlineJavascriptRequirement
- class: StepInputExpressionRequirement
- class: MultipleInputFeatureRequirement

inputs:
- id: inp1
  type:
  - string
  - 'null'
- id: inp2
  type:
  - string
  - 'null'

outputs:
- id: out
  type:
    type: array
    items: File
  outputSource: print/outs

steps:
- id: print
  in:
  - id: _print_inps_inp1
    source: inp1
  - id: _print_inps_inp2
    source: inp2
  - id: inps
    valueFrom: |-
      $([(inputs._print_inps_inp1 != null) ? inputs._print_inps_inp1 : "default1", (inputs._print_inps_inp2 != null) ? (inputs._print_inps_inp2 + "_suffix") : ""])
  run: tools/ArrayStepTool.cwl
  out:
  - id: outs
id: cwl_test_array_step_input
"""
