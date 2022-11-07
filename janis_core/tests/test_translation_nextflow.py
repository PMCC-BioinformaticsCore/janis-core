from copy import deepcopy
import unittest
from typing import Optional
import regex as re

from janis_core.tests.testtools import (
    InputQualityTestTool,
    EchoTestTool,
    FilenameGeneratedTool,
    OperatorResourcesTestTool,
    SplitTextTestTool,
    SumTestTool,
    JoinArrayTestTool,
    FileInputTestTool,
    SecondaryInputTestTool,
    BasicTestTool,
    AppendedSecondaryOutputTestTool,
    ReplacedSecondaryOutputTestTool,
)

from janis_core.tests.testworkflows import (
    StepInputExpressionTestWF,
    ConditionStepTestWF,
    ArraysOfSecondaryFilesOutputsTestWF,
    AliasSelectorTestWF,
)

from janis_core import (
    WorkflowBuilder,
    ToolInput,
    InputSelector,
    StringFormatter,
    JoinOperator,
)

from janis_core.translations import NextflowTranslator as translator
from janis_core.translations import nfgen
from janis_core.translations.nfgen import settings
from janis_core import Array, String, File, Boolean, Filename



class TestSettings(unittest.TestCase):

    def test_get_settings(self):
        self.assertEquals(settings.LIB_FILENAME, 'lib.nf')
        self.assertEquals(settings.CONFIG_FILENAME, 'nextflow.config')
    
    def test_set_settings(self):
        settings.MINIMAL_PROCESS = False
        self.assertEquals(settings.MINIMAL_PROCESS, False)

    @unittest.skip('not implemented')
    def test_minimal_process(self):
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_janis_assistant_process(self):
        raise NotImplementedError
    


WF_INPUTS_SINGLES = {
    'testinput': 'null',
    'fastqc1_adapters': 'null',
    'fastqc1_contaminants': 'null',
    'fastqc1_limits': 'null',
    'fastqc2_adapters': 'null',
    'fastqc2_contaminants': 'null',
    'fastqc2_limits': 'null',
    'prokka_proteins': 'null',
}

WF_INPUTS_ARRS = {
    'inforwardreads': 'null',
    'inlongreads': 'null',
    'inreversereads': 'null',
}

TOOL_INPUTS = {
    'prokka_genus': '"Escherichia"',
    'prokka_increment': '10',
    'prokka_locustag': '"PROKKA"',
    'prokka_species': '"Coli"',
    'prokka_strain': '"C-1"',
    'quast_contigthresholds': 'null',
    'quast_referenceslist': '"temp_ref_list_fp"',
    'unicycler_kmers': '""',
    'unicycler_scores': '""',
    'unicycler_startgenecov': '95.0',
    'unicycler_startgeneid': '90.0',
}

COMBINED_INPUTS = WF_INPUTS_SINGLES | WF_INPUTS_ARRS | TOOL_INPUTS

class TestParamRegistration(unittest.TestCase):
    MINIMAL_PARAMS = {
        'inforwardreads': 'null',
        'inreversereads': 'null',
        'inlongreads': 'null',
        'testinput': 'null',
        'fastqc1_adapters': 'null',
        'fastqc1_contaminants': 'null',
        'fastqc1_limits': 'null',
        'fastqc2_adapters': 'null',
        'fastqc2_contaminants': 'null',
        'fastqc2_limits': 'null',
        'unicycler_kmers': '""',
        'unicycler_scores': '""',
        'unicycler_startgenecov': '95.0',
        'unicycler_startgeneid': '90.0',
    }
    
    ADDITIONAL_PARAMS = {
        'fastqc1_extract': 'null',
        'fastqc1_nogroup': 'null',
        'fastqc1_quiet': 'null',
        'fastqc1_kmers': 'null',
        'fastqc1_minlength': 'null',
        'fastqc1_optionf': 'null',
        'fastqc1_outdir': 'null',
        'fastqc2_extract': 'null',
        'fastqc2_nogroup': 'null',
        'fastqc2_quiet': 'null',
        'fastqc2_kmers': 'null',
        'fastqc2_minlength': 'null',
        'fastqc2_optionf': 'null',
        'fastqc2_outdir': 'null',
        'fastqc3_adapters': 'null',
        'fastqc3_contaminants': 'null',
        'fastqc3_limits': 'null',
        'fastqc3_extract': 'null',
        'fastqc3_nogroup': 'null',
        'fastqc3_quiet': 'null',
        'fastqc3_kmers': 'null',
        'fastqc3_minlength': 'null',
        'fastqc3_optionf': 'null',
        'fastqc3_outdir': 'null',
        'unicycler_largestcomponent': 'null',
        'unicycler_nocorrect': 'null',
        'unicycler_nopilon': 'null',
        'unicycler_norotate': 'null',
        'unicycler_contamination': 'null',
        'unicycler_depthfilter': 'null',
        'unicycler_kmercount': 'null',
        'unicycler_linearseqs': 'null',
        'unicycler_lowscore': 'null',
        'unicycler_maxkmerfrac': 'null',
        'unicycler_minanchorseglen': 'null',
        'unicycler_mincomponentsize': 'null',
        'unicycler_mindeadendsize': 'null',
        'unicycler_minfastalength': 'null',
        'unicycler_minkmerfrac': 'null',
        'unicycler_minpolishsize': 'null',
        'unicycler_mode': 'null',
        'unicycler_optiono': 'null',
        'unicycler_options': 'null',
        'unicycler_optiont': 'null',
        'unicycler_pilonpath': 'null',
        'unicycler_startgenes': 'null',
        'unicycler_verbosity': 'null',
    }

    FULL_PARAMS = MINIMAL_PARAMS | ADDITIONAL_PARAMS

    UNICYCLER_MINIMAL = {
        'kmers': '""',
        'scores': '""',
        'startgenecov': '95.0',
        'startgeneid': '90.0',
    }

    UNICYCLER_ADDITIONAL = {
        'largestcomponent': 'null',
        'nocorrect': 'null',
        'nopilon': 'null',
        'norotate': 'null',
        'contamination': 'null',
        'depthfilter': 'null',
        'kmercount': 'null',
        'linearseqs': 'null',
        'lowscore': 'null',
        'maxkmerfrac': 'null',
        'minanchorseglen': 'null',
        'mincomponentsize': 'null',
        'mindeadendsize': 'null',
        'minfastalength': 'null',
        'minkmerfrac': 'null',
        'minpolishsize': 'null',
        'mode': 'null',
        'optiono': 'null',
        'options': 'null',
        'optiont': 'null',
        'pilonpath': 'null',
        'startgenes': 'null',
        'verbosity': 'null',
    }

    UNICYCLER_FULL = UNICYCLER_MINIMAL | UNICYCLER_ADDITIONAL

    def setUp(self) -> None:
        from janis_core.tests.testworkflows.assembly import w
        from janis_core.tests.data.janis.simple_truncated.tools.unicycler import unicycler
        self.wf = w
        self.unicycler = unicycler
        self.maxDiff = None

    def refresh_params_workflowmode(self) -> None:
        scope: list[str] = []
        nfgen.params.clear()
        nfgen.params.register(the_entity=self.wf, scope=scope)
        for step in self.wf.step_nodes.values():
            current_scope = deepcopy(scope)
            current_scope.append(step.id())
            nfgen.params.register(the_entity=step.tool, sources=step.sources, scope=current_scope)
    
    def refresh_params_toolmode(self) -> None:
        nfgen.params.clear()
        nfgen.params.register(the_entity=self.unicycler, sources={}, scope=[])

    def test_all_full_workflowmode(self) -> None:
        # check each necessary param is registered for 
        # full process translation in workflow mode
        settings.MINIMAL_PROCESS = False
        settings.MODE = 'workflow'
        self.refresh_params_workflowmode()
        params = nfgen.params.getall()
        params_dict = {p.name: p.value for p in params}
        self.assertEquals(params_dict, self.FULL_PARAMS)
    
    def test_all_minimal_workflowmode(self) -> None:
        # check each necessary param is registered for 
        # minimal process translation in workflow mode
        settings.MINIMAL_PROCESS = True
        settings.MODE = 'workflow'
        self.refresh_params_workflowmode()
        params = nfgen.params.getall()
        params_dict = {p.name: p.value for p in params}
        self.assertEquals(params_dict, self.MINIMAL_PARAMS)
    
    @unittest.skip('not implemented')
    def test_all_full_toolmode(self) -> None:
        # check each necessary param is registered for 
        # full process translation in tool mode
        settings.MINIMAL_PROCESS = False
        settings.MODE = 'tool'
        self.refresh_params_toolmode()
        params = nfgen.params.getall()
        params_dict = {p.name: p.value for p in params}
        self.assertEquals(params_dict, self.MINIMAL_PARAMS)
    
    def test_all_minimal_toolmode(self) -> None:
        # check each necessary param is registered for 
        # minimal process translation in tool mode
        settings.MINIMAL_PROCESS = True
        settings.MODE = 'tool'
        self.refresh_params_toolmode()
        params = nfgen.params.getall()
        params_dict = {p.name: p.value for p in params}
        self.assertEquals(params_dict, self.UNICYCLER_MINIMAL)

    @unittest.skip('not implemented')
    def test_workflow_inputs(self) -> None:
        raise NotImplementedError
        # check each wf input now has a param
        params = nfgen.params.getall()
            
        wf_singles = {}
        for p in params:
            if p.is_wf_input and not isinstance(p.dtype, Array):
                wf_singles[p.name] = p.value

        self.assertEquals(wf_singles, WF_INPUTS_SINGLES)
    
    @unittest.skip('not implemented')
    def test_tool_inputs(self) -> None:
        raise NotImplementedError
        params = nfgen.params.getall()
            
        exposed_inputs = {}
        for p in params:
            if not p.is_wf_input:
                exposed_inputs[p.name] = p.value

        self.assertEquals(exposed_inputs, TOOL_INPUTS)



class TestNextflowConfig(unittest.TestCase):

    def setUp(self) -> None:
        from janis_core.tests.data.janis.simple.workflow import w
        scope: list[str] = []
        nfgen.params.clear()
        nfgen.params.register(the_entity=w, scope=scope)
        for step in w.step_nodes.values():
            current_scope = deepcopy(scope)
            current_scope.append(step.id())
            nfgen.params.register(the_entity=step.tool, sources=step.sources, scope=current_scope)

    def test_workflow_config(self) -> None:
        config = translator.stringify_translated_inputs({})
        # basic structure
        self.assertIn('docker.enabled = true', config)
        self.assertIn('params {\n\n', config)
        self.assertIn('\n\n}', config)
        # expected params are present & correct
        for name, val in COMBINED_INPUTS.items():
            pattern = f'{name}.*?{val}'
            matches = re.findall(pattern, config)
            self.assertGreater(len(matches), 0)

    @unittest.skip('not implemented')
    def test_tool_config(self) -> None:
        raise NotImplementedError
    


class TestWorkflowInputs(unittest.TestCase):

    def setUp(self) -> None:
        from janis_core.tests.data.janis.simple.workflow import w
        self.wf = w
        self.inputs = list(self.wf.input_nodes.values())
        self.file_inputs = set([inp for inp in self.inputs if isinstance(inp.datatype, File)])
        self.channels = translator.gen_channel_declarations(self.wf).channels
    
    def test_params(self) -> None:
        """
        Every wf input should have a param. 
        """
        params = translator.build_inputs_file(self.wf)
        for winp in self.inputs:
            self.assertIn(winp.id(), params)
    
    def test_file_channels(self) -> None:
        """
        Every File type wf input should have a channel. 
        """
        channel_src_names = set([c.wfinp_name for c in self.channels])
        for winp in self.file_inputs:
            self.assertIn(winp.id(), channel_src_names)
    
    def test_optional_file_channels(self) -> None:
        """
        Every Optional(File) type wf input should have a channel. 
        '.ifEmpty(null)' should appear in the channel string definition.
        """
        optional_inputs = [x for x in self.file_inputs if x.datatype.optional == True]
        for winp in optional_inputs:
            channel = [c for c in self.channels if c.wfinp_name == winp.id()][0]
            self.assertIn('.ifEmpty(null)', channel.get_string())

    def test_nonfile_no_channel(self) -> None:
        """
        Non-File-type wf input should not have channels.
        """
        nonfile_inputs = [x for x in self.inputs if not isinstance(x.datatype, File)]
        channel_src_names = set([c.wfinp_name for c in self.channels])
        for winp in nonfile_inputs:
            self.assertNotIn(winp.id(), channel_src_names)
    



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

        workflow = ConditionStepTestWF()
        step_keys = list(workflow.step_nodes.keys())

        step_id = "print"
        tool = workflow.step_nodes[step_id].tool
        inputs = translator.gen_step_inval_dict(tool)
        expected = {"inp": "[$params.mystring, $get_string.out.out].first()"}

        self.assertEqual(expected, inputs)

    def test_simple(self):
        w1 = ArraysOfSecondaryFilesOutputsTestWF()
        w1_step_keys = list(w1.step_nodes.keys())

        expected = {"testtool": "$params.inp"}
        self.assertEqual(
            expected,
            translator.gen_step_inval_dict(w1.step_nodes["stp"].tool),
        )

    def test_with_expression(self):
        w2 = StepInputExpressionTestWF()
        w2_step_keys = list(w2.step_nodes.keys())

        expected = {
            "inp": "$params.mystring ? $params.mystring : $params.mystring_backup"
        }
        self.assertEqual(
            expected,
            translator.gen_step_inval_dict(w2.step_nodes["print"].tool),
        )

    def test_multi_steps(self):
        w3 = AliasSelectorTestWF()
        w3_step_keys = list(w3.step_nodes.keys())

        expected1 = {"testtool": "$params.inp"}
        self.assertEqual(
            expected1,
            translator.gen_step_inval_dict(w3.step_nodes["stp1"].tool),
        )

        expected2 = {"inp": "$stp1.out.out"}
        self.assertEqual(
            expected2,
            translator.gen_step_inval_dict(w3.step_nodes["stp2"].tool),
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
    w1 = ArraysOfSecondaryFilesOutputsTestWF()
    w2 = StepInputExpressionTestWF()
    w3 = AliasSelectorTestWF()

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
    any_tool = BasicTestTool()

    def test_string_formatter(self):
        b = StringFormatter("no format")
        res = nfgen.translate_string_formatter(b, self.any_tool)
        self.assertEqual("no format", res)

    def test_string_formatter_one_string_param(self):
        b = StringFormatter("there's {one} arg", one="a string")
        res = nfgen.translate_string_formatter(b, self.any_tool)
        self.assertEqual("there's ${'a string'} arg", res)

    def test_string_formatter_one_input_selector_param(self):
        b = StringFormatter("an input {arg}", arg=InputSelector("testtool"))
        res = nfgen.translate_string_formatter(
            b, self.any_tool, inputs_dict=self.any_tool.inputs_map()
        )
        self.assertEqual("an input ${testtool}", res)

    def test_string_formatter_two_param(self):
        tool = InputQualityTestTool()
        b = StringFormatter(
            "{username}:{password}",
            username=InputSelector("user"),
            password=InputSelector("static"),
        )
        res = nfgen.translate_string_formatter(
            b, tool, inputs_dict=tool.inputs_map()
        )
        self.assertEqual(
            "${user}:${static}",
            res,
        )

    def test_escaped_characters(self):
        tool = InputQualityTestTool()
        b = StringFormatter(
            "{username}\\t{password}",
            username=InputSelector("user"),
            password=InputSelector("static"),
        )
        res = nfgen.translate_string_formatter(
            b, tool, inputs_dict=tool.inputs_map()
        )
        self.assertEqual("${user}\\t${static}", res)

        res2 = nfgen.translate_string_formatter(
            b, tool, inputs_dict=tool.inputs_map(), in_shell_script=True
        )
        self.assertEqual("${user}\\\\t${static}", res2)

    def test_expression_arg(self):
        tool = BasicTestTool()
        b = StringFormatter(
            "{name}:{items}",
            name=InputSelector("testtool"),
            items=JoinOperator(InputSelector("arrayInp"), separator=";"),
        )

        res = nfgen.translate_string_formatter(
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
            {"inp": f"/{settings.NO_FILE_PATH_PREFIX}1"},
            self.translator.build_inputs_file(wf),
        )

    def test_file_type_no_secondary_no_value(self):
        wf = WorkflowBuilder("test_file_input")
        wf.input("inp", DataTypeNoSecondary())
        self.assertDictEqual(
            {"inp": f"/{settings.NO_FILE_PATH_PREFIX}1"},
            self.translator.build_inputs_file(wf),
        )

    def test_file_type_no_secondary_multiple_inputs(self):
        wf = WorkflowBuilder("test_file_input")
        wf.input("inp", DataTypeNoSecondary())
        wf.input("inp2", DataTypeNoSecondary())
        self.assertDictEqual(
            {
                "inp": f"/{settings.NO_FILE_PATH_PREFIX}1",
                "inp2": f"/{settings.NO_FILE_PATH_PREFIX}2",
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
    settings.MODE = 'tool'
    settings.MINIMAL_PROCESS = False
    maxDiff: Optional[int] = None
    
    def test_stdout_out_tool(self):
        scope = []
        values = {}
        tool = EchoTestTool()
        p = translator.gen_process_from_cmdtool(tool, values, scope)
        expected_contents = [
            'process EchoTestTool',
            'publishDir "$params.outdir"',
            'output:',
            'path "janisstdout_EchoTestTool" , emit: out',
            'script:',
            'echo',
            '${params.inp}'
        ]
        for line in expected_contents:
            self.assertIn(line, p.get_string())

    def test_operator_resource_tool(self):
        scope = []
        values = {}
        tool = OperatorResourcesTestTool()
        p = translator.gen_process_from_cmdtool(tool, values, scope)
        # f"""
        # process EchoTestTool
        # {{
        #   input:
        #     path inputFile
        #     val outputFiles

        #   output:
        #     path "${{'janisstdout_EchoTestTool'}}" , emit: out

        #   publishDir "$launchDir/EchoTestTool"
        #   memory "${{params.runtime_memory ? params.runtime_memory + 'GB': ''}}"
        #   cpus "${{params.runtime_cpu ? params.runtime_cpu : ''}}"
        #   disk "${{params.runtime_disks ? params.runtime_disks + 'GB': ''}}"
        #   time "${{params.runtime_seconds + 's'}}"

        #   script:

        #     def inputFileWithPrefix = apply_prefix(inputFile, ' ', 'False')

        #     def outputFilesWithPrefix = apply_prefix(outputFiles, ' ', 'False')

        #     def runtime_memory = params.runtime_memory

        #     def runtime_cpu = params.runtime_cpu

        #     def runtime_disks = params.runtime_disks

        #     def runtime_seconds = params.runtime_seconds
        #     \"\"\"
        #     echo \\
        #     $inputFileWithPrefix | tee janisstdout_EchoTestTool
        #     \"\"\"
        # }}
        # """
        expected_contents = [
            'process EchoTestTool',
            'input:',
            'output:',
        ]
        # TODO ask richard - unnecessary?
        raise NotImplementedError
        print(p.get_string())
        for line in expected_contents:
            self.assertIn(line, p.get_string())
        
    def test_filename_generated_tool(self):
        # TODO
        scope = []
        values = {}
        tool = FilenameGeneratedTool()
        p = translator.gen_process_from_cmdtool(tool, values, scope)
        print(p.get_string())

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
        scope = []
        values = {}
        tool = SecondaryInputTestTool()
        p = translator.gen_process_from_cmdtool(tool, values, scope)
        print(p.get_string())
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
        scope = []
        values = {}
        tool = AppendedSecondaryOutputTestTool()
        p = translator.gen_process_from_cmdtool(tool, values, scope)
        print(p.get_string())
        expected_contents = [
            'process TestTranslationtool',
            'tuple path("${testtool + ".bam"}"), path("${testtool + ".bam.bai"}") , emit: out',
            'def arrayInp = params.arrayInp ? "${params.arrayInp}" : ""',
            'echo',
            '${params.testtool}',
            '${arrayInp}',
            'test:\\t:escaped:\\n:characters"',
        ]
        for line in expected_contents:
            self.assertIn(line, p.get_string())


    def test_tool_with_replaced_secondary_output(self):
        p = translator.gen_process_from_cmdtool(
            ReplacedSecondaryOutputTestTool()
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
        p = translator.gen_process_from_codetool(SplitTextTestTool())
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
        p = translator.gen_process_from_codetool(SumTestTool())
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
        p = translator.gen_process_from_codetool(JoinArrayTestTool())
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
        p = translator.gen_process_from_codetool(FileInputTestTool())
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
            SecondaryInputTestTool()
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
