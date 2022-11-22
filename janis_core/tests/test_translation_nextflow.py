from copy import deepcopy
import unittest
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
    CatTestTool,
    SecondaryInputTestTool,
    BasicTestTool,
    AppendedSecondaryOutputTestTool,
    ReplacedSecondaryOutputTestTool,
)

from janis_core.tests.testworkflows import (

    # basics
    BasicIOTestWF,
    WildcardSelectorOutputTestWF,
    InputSelectorOutputTestWF,
    StepInputsTestWF,
    StepInputsWFInputTestWF,
    StepInputsStaticTestWF,
    StepInputsPartialStaticTestWF,
    StepInputsMinimalTestWF,
    StepConnectionsTestWF,
    DirectivesTestWF,

    # arrays
    ArrayIOTestWF,
    ArrayIOExtrasTestWF,
    ArrayStepInputsTestWF,
    ArrayStepConnectionsTestWF,

    # scatter
    BasicScatterTestWF,
    ChainedScatterTestWF,
    ScatterDotTestWF,
    ScatterCrossTestWF,

    # secondaries
    SecondariesIOTestWF,
    SecondariesConnectionsTestWF,

    # combos
    ScatterSecondariesTestWF,
    ArraySecondariesTestWF,

    # additional features
    StepInputExpressionTestWF,
    ConditionStepTestWF,
    AliasSelectorTestWF,
    ArraysOfSecondaryFilesOutputsTestWF,
    ForEachTestWF,

    # specific workflows
    AssemblyTestWF,
)

from janis_core import (
    WorkflowBuilder,
    Workflow,
    CommandTool,
    InputSelector,
    StringFormatter,
    JoinOperator,
)

from janis_core.translations import NextflowTranslator as translator
from janis_core.translations import nfgen
from janis_core.translations.nfgen import settings
from janis_core.translations.nfgen.params import Param
from janis_core.translations.nfgen.channels import Channel
from janis_core import (
    Array, 
    String, 
    Int, 
    Float, 
    File, 
    Boolean, 
    Filename
)
from janis_core.translations.nfgen.utils import to_groovy_str


### helper classes


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



### helper functions


def refresh_workflow_inputs(wf: Workflow) -> None:
    nfgen.params.clear()
    nfgen.channels.clear()
    nfgen.register_workflow_inputs(wf, scope=[])

# def refresh_workflow_inputs(wf: Workflow) -> dict[str, str]:
#     scope: list[str] = []
#     nfgen.params.clear()
#     nfgen.params.register(the_entity=wf, scope=scope)
#     # for step in wf.step_nodes.values():
#     #     current_scope = deepcopy(scope)
#     #     current_scope.append(step.id())
#     #     nfgen.params.register(the_entity=step.tool, sources=step.sources, scope=current_scope)
#     params = nfgen.params.getall()
#     return {p.name: p.groovy_value for p in params} 

# def register_params_tool(tool: CommandTool) -> dict[str, str]:
#     scope: list[str] = []
#     nfgen.params.clear()
#     nfgen.params.register(the_entity=tool, scope=scope)
#     params = nfgen.params.getall()
#     return {p.name: p.groovy_value for p in params} 

# def register_channels(wf: Workflow) -> list[Channel]:
#     nfgen.channels.clear()
#     nfgen.channels.register(wf)
#     return nfgen.channels.getall()



class TestToGroovyStr(unittest.TestCase):

    def test_string(self) -> None:
        inp = 'Hello'
        expected = "'Hello'"
        actual = to_groovy_str(inp, String())
        self.assertEquals(expected, actual)
    
    def test_numeric(self) -> None:
        # int
        inp = 5
        expected = '5'
        actual = to_groovy_str(inp, Int())
        self.assertEquals(expected, actual)
        # float
        inp = 5.2
        expected = '5.2'
        actual = to_groovy_str(inp, Float())
        self.assertEquals(expected, actual)
    
    def test_bool(self) -> None:
        # true
        inp = True
        expected = 'true'
        actual = to_groovy_str(inp, Boolean())
        self.assertEquals(expected, actual)
        # false
        inp = False
        expected = 'false'
        actual = to_groovy_str(inp, Boolean())
        self.assertEquals(expected, actual)
    
    def test_none(self) -> None:
        inp = None
        expected = 'null'
        actual = to_groovy_str(inp, String())
        self.assertEquals(expected, actual)
    
    def test_empty_array(self) -> None:
        inp = []
        expected = "[]"
        actual = to_groovy_str(inp, String())
        self.assertEquals(expected, actual)
    
    def test_string_array(self) -> None:
        inp = ['Hello']
        expected = "['Hello']"
        actual = to_groovy_str(inp, String())
        self.assertEquals(expected, actual)
    
    def test_numeric_array(self) -> None:
        # int
        inp = [1, 2, 3]
        expected = '[1, 2, 3]'
        actual = to_groovy_str(inp, Int())
        self.assertEquals(expected, actual)
        # float
        inp = [1.1, 2.2, 3.3]
        expected = '[1.1, 2.2, 3.3]'
        actual = to_groovy_str(inp, Float())
        self.assertEquals(expected, actual)
    
    def test_bool_array(self) -> None:
        inp = [True, False]
        expected = '[true, false]'
        actual = to_groovy_str(inp)
        self.assertEquals(expected, actual)
    
    def test_none_array(self) -> None:
        inp = [None, None]
        expected = '[null, null]'
        actual = to_groovy_str(inp)
        self.assertEquals(expected, actual)




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
    



class TestParams(unittest.TestCase):
    
    def setUp(self) -> None:
        settings.MINIMAL_PROCESS = True
        settings.MODE = 'workflow'
    
    def test_channel_params(self) -> None:
        """
        Every channel requires a param. 
        Channels are created from a subset of workflow inputs.
        """
        wf = AssemblyTestWF()
        refresh_workflow_inputs(wf)
        params = nfgen.params.getall()
        param_ids = {p.name for p in params}
        expected_ids = {
            'unicycler_start_gene_id', 
            'test_input', 
            'fastqc2_limits', 
            'unicycler_start_gene_cov', 
            'fastqc1_adapters', 
            'unicycler_scores', 
            'unicycler_kmers', 
            'fastqc2_adapters', 
            'fastqc2_contaminants', 
            'in_long_reads', 
            'in_reverse_reads', 
            'in_forward_reads', 
            'fastqc1_contaminants', 
            'fastqc1_limits'
        }
        for inp_id in expected_ids:
            self.assertIn(inp_id, param_ids)

    def test_fed_step_inputs(self) -> None:
        wf = StepInputsWFInputTestWF()
        
        # full
        # settings.MINIMAL_PROCESS = False
        # actual_params = refresh_workflow_inputs(wf)
        # expected_params = {
        #     'in_file': 'null',
        #     'stp1_pos_default': 'null',
        #     'stp1_pos_optional': 'null',
        #     'stp1_flag_true': 'null',
        #     'stp1_flag_false': 'null',
        #     'stp1_opt_basic': 'null',
        #     'stp1_opt_default': 'null',
        #     'stp1_opt_optional': 'null'
        # }
        # self.assertEquals(actual_params, expected_params)
        
        # minimal
        settings.MINIMAL_PROCESS = True
        refresh_workflow_inputs(wf)
        actual_params = nfgen.params.serialize()
        expected_params = {
            'in_file': 'null',
            'in_file_opt': 'null',
            'in_str': 'null',
            'in_int': 'null',
            'in_bool': 'null',
        }
        self.assertEquals(actual_params, expected_params)

    def test_static_step_inputs(self) -> None:
        wf = StepInputsStaticTestWF()
        
        # full
        # settings.MINIMAL_PROCESS = False
        # actual_params = refresh_workflow_inputs(wf)
        # expected_params = {
        #     'in_file': 'null',
        #     'stp2_pos_default': '100',
        #     'stp2_pos_optional': "'static'",
        #     'stp2_flag_true': 'false',
        #     'stp2_flag_false': 'true',
        #     'stp2_opt_basic': "'static'",
        #     'stp2_opt_default': '100',
        #     'stp2_opt_optional': "''"
        # }
        # self.assertEquals(actual_params, expected_params)
        
        # minimal
        settings.MINIMAL_PROCESS = True
        refresh_workflow_inputs(wf)
        actual_params = nfgen.params.serialize()
        expected_params = {
            'in_file': 'null',
            'in_str': 'null',
            'in_int': 'null',
            'in_bool': 'null',
            'stp2_pos_default': '100',
            'stp2_pos_optional': "'static'",
            'stp2_flag_true': 'false',
            'stp2_flag_false': 'true',
            'stp2_opt_basic': "'static'",
            'stp2_opt_default': '100',
            'stp2_opt_optional': "''"
        }
        self.assertEquals(actual_params, expected_params)
    
    def test_static_and_omitted_step_inputs(self) -> None:
        wf = StepInputsPartialStaticTestWF()
        
        # # full
        # settings.MINIMAL_PROCESS = False
        # actual_params = refresh_workflow_inputs(wf)
        # expected_params = {
        #     'in_file': 'null',
        #     'stp3_pos_default': 'null',
        #     'stp3_pos_optional': 'null',
        #     'stp3_flag_true': 'null',
        #     'stp3_flag_false': 'null',
        #     'stp3_opt_basic': "'static'",
        #     'stp3_opt_default': '100',
        #     'stp3_opt_optional': "''"
        # }
        # self.assertEquals(actual_params, expected_params)
        
        # minimal
        settings.MINIMAL_PROCESS = True
        refresh_workflow_inputs(wf)
        actual_params = nfgen.params.serialize()
        expected_params = {
            'in_file': 'null',
            'in_str': 'null',
            'in_int': 'null',
            'in_bool': 'null',
            'stp3_opt_basic': "'static'",
            'stp3_opt_default': '100',
            'stp3_opt_optional': "''"
        }
        self.assertEquals(actual_params, expected_params)
    
    def test_omitted_step_inputs(self) -> None:
        wf = StepInputsMinimalTestWF()
        
        # # full
        # settings.MINIMAL_PROCESS = False
        # actual_params = refresh_workflow_inputs(wf)
        # expected_params = {
        #     'in_file': 'null',
        #     'stp4_pos_default': 'null',
        #     'stp4_pos_optional': 'null',
        #     'stp4_flag_true': 'null',
        #     'stp4_flag_false': 'null',
        #     'stp4_opt_basic': 'null',
        #     'stp4_opt_default': 'null',
        #     'stp4_opt_optional': 'null'
        # }
        # self.assertEquals(actual_params, expected_params)
        
        # minimal
        settings.MINIMAL_PROCESS = True
        refresh_workflow_inputs(wf)
        actual_params = nfgen.params.serialize()
        expected_params = {
            'in_file': 'null',
            'in_str': 'null',
            'in_int': 'null',
            'in_bool': 'null',
        }
        self.assertEquals(actual_params, expected_params)

    def test_array_inputs(self) -> None:
        wf = ArrayIOTestWF()
        refresh_workflow_inputs(wf)
        actual_params = nfgen.params.serialize()
        expected_params = {
            'in_file_array': '[]',
            'in_str_array': '[]',
            'in_int_array': '[]',
        }
        self.assertEquals(actual_params, expected_params)

    def test_secondaries_inputs(self) -> None:
        wf = SecondariesIOTestWF()
        refresh_workflow_inputs(wf)
        actual_params = nfgen.params.serialize()
        expected_params = {
            'in_alignments_bam': 'null',
            'in_alignments_bai': 'null',
        }
        self.assertEquals(actual_params, expected_params)
    
    def test_secondaries_array_inputs(self) -> None:
        wf = ArraySecondariesTestWF()
        refresh_workflow_inputs(wf)
        actual_params = nfgen.params.serialize()
        expected_params = {
            'in_alignments_bams': '[]',
            'in_alignments_bais': '[]',
        }
        self.assertEquals(actual_params, expected_params)
    
    @unittest.skip('not implemented')
    def test_translate_commandtool(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_translate_codetool(self) -> None:
        raise NotImplementedError




class TestChannels(unittest.TestCase):
    
    def setUp(self) -> None:
        settings.MINIMAL_PROCESS = True
        settings.MODE = 'workflow'
        self.maxDiff = None

    def test_infer_wf_inputs(self) -> None:
        # checks the inferred wf inputs (from total wf inputs) are correct
        wf = AssemblyTestWF()
        refresh_workflow_inputs(wf)
        channel_ids = {c.name for c in nfgen.channels.getall()}
        expected_ids = {
            'ch_in_forward_reads',
            'ch_in_reverse_reads',
            'ch_in_long_reads',
            'ch_test_input',
            'ch_fastqc1_adapters',
            'ch_fastqc1_contaminants',
            'ch_fastqc1_limits',
            'ch_fastqc2_adapters',
            'ch_fastqc2_contaminants',
            'ch_fastqc2_limits',
        }
        self.assertEqual(channel_ids, expected_ids)
    
    def test_optional_file_channels(self) -> None:
        """
        Every Optional(File) type wf input should have a channel. 
        '.ifEmpty(null)' should appear in the channel string definition.
        """
        wf = AssemblyTestWF()
        refresh_workflow_inputs(wf)
        relevant_channel_names = {
            'ch_fastqc1_adapters',
            'ch_fastqc1_contaminants',
            'ch_fastqc1_limits',
            'ch_fastqc2_adapters',
            'ch_fastqc2_contaminants',
            'ch_fastqc2_limits',
        }
        for name in relevant_channel_names:
            channel = nfgen.channels.get(name)
            self.assertTrue(channel)
            self.assertIn('.ifEmpty( null )', channel.get_string())

    def test_nonfile_no_channel(self) -> None:
        """
        Non-File-type wf input should not have channels.
        """
        wf = AssemblyTestWF()
        refresh_workflow_inputs(wf)
        nonfile_wf_input_ids = {
            'unicycler_kmers',
            'unicycler_scores',
            'unicycler_startGeneCov',
            'unicycler_startGeneId',
        }
        channel_janis_references = {c.ref_name for c in nfgen.channels.getall()}
        for winp_id in nonfile_wf_input_ids:
            self.assertNotIn(winp_id, channel_janis_references)

    def test_array_inputs(self) -> None:
        wf = ArrayIOTestWF()
        refresh_workflow_inputs(wf)
        # channels created
        channels_ids = {c.name for c in nfgen.channels.getall()}
        expected_ids = {
            'ch_in_file_array',
        }
        self.assertEqual(channels_ids, expected_ids)
        # channels are marked to collect
        for c in nfgen.channels.getall():
            self.assertTrue(c.collect)
        
    def test_secondaries_inputs(self) -> None:
        wf = SecondariesIOTestWF()
        refresh_workflow_inputs(wf)
        channels_ids = {c.name for c in nfgen.channels.getall()}
        expected_ids = {
            'ch_in_alignments',
        }
        self.assertEqual(channels_ids, expected_ids)
        alignments_ch = nfgen.channels.get('ch_in_alignments')
        self.assertTrue(len(alignments_ch.params) == 2)
        self.assertTrue(alignments_ch.collect)
    
    def test_secondaries_array_inputs(self) -> None:
        wf = ArraySecondariesTestWF()
        refresh_workflow_inputs(wf)
        channels_ids = {c.name for c in nfgen.channels.getall()}
        expected_ids = {
            'ch_in_alignments_bams',
            'ch_in_alignments_bais',
        }
        self.assertEqual(channels_ids, expected_ids)
        for channel in nfgen.channels.getall():
            self.assertTrue(len(channel.params) == 1)
            self.assertTrue(channel.collect)

    @unittest.skip('not implemented')
    def test_channel_methods(self) -> None:
        raise NotImplementedError

    @unittest.skip('not implemented')
    def test_translate_commandtool(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_translate_codetool(self) -> None:
        raise NotImplementedError




class TestProcessDirectives(unittest.TestCase):
    """
    Tests identifying tool inputs which should be process inputs.
    Need a process input for each tool input in step sources

    INCLUDES WORKFLOW OUTPUTS as publishDir
    """

    def setUp(self) -> None:
        settings.MINIMAL_PROCESS = True
        settings.MODE = 'workflow'

    def test_directives_order(self) -> None:
        wf = DirectivesTestWF()
        refresh_workflow_inputs(wf)
        tool = wf.step_nodes["stp1"].tool
        sources = wf.step_nodes["stp1"].sources
        process = translator.gen_process_from_cmdtool(tool, sources, scope=['stp1'])
        process = translator.handle_container(tool, process)
        directives = nfgen.process.order_directives(process.directives)
        actual_order = [type(x).__name__ for x in directives]
        expected_order = [
            'ContainerDirective',
            'PublishDirDirective',
            'DebugDirective',
            'DiskDirective',
            'MemoryDirective',
            'TimeDirective',
        ]
        for actual, expected in zip(actual_order, expected_order):
            self.assertEqual(actual, expected)

    def test_directives(self) -> None:
        wf = DirectivesTestWF()
        refresh_workflow_inputs(wf)
        tool = wf.step_nodes["stp1"].tool
        sources = wf.step_nodes["stp1"].sources
        process = translator.gen_process_from_cmdtool(tool, sources, scope=['stp1'])
        process = translator.handle_container(tool, process)
        actual_directives = {d.get_string() for d in process.directives}
        expected_directives = {
            'container "quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0"',
            'publishDir "${params.outdir}/stp1"',
            'debug true',
            'disk "${params.stp1_disk}"',
            'memory "${params.stp1_memory}"',
            'time "${params.stp1_time}"'
        }
        for direc in expected_directives:
            self.assertIn(direc, actual_directives)
    
    @unittest.skip('not implemented')
    def test_translate_commandtool(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_translate_codetool(self) -> None:
        raise NotImplementedError
    



class TestProcessInputs(unittest.TestCase):
    """
    Tests identifying tool inputs which should be process inputs.
    Need a process input for each tool input in step sources
    """

    def setUp(self) -> None:
        settings.MINIMAL_PROCESS = True
        settings.MODE = 'workflow'

    def test_wf_inputs(self) -> None:
        # need a process input for each File wf input in step sources.
        # non-files are fed data via params.
        wf = AssemblyTestWF()
        refresh_workflow_inputs(wf)
        tool = wf.step_nodes["unicycler"].tool
        sources = wf.step_nodes["unicycler"].sources
        process = translator.gen_process_from_cmdtool(tool, sources, scope=['unicycler'])
        expected_inputs = {
            'path option1',
            'path option2',
            'path optionL',
        }
        actual_inputs = {inp.get_string() for inp in process.inputs}
        self.assertEqual(actual_inputs, expected_inputs)
    
    def test_connections(self) -> None:
        # need a process input for each connection in step sources.
        wf = AssemblyTestWF()
        refresh_workflow_inputs(wf)
        tool = wf.step_nodes["CatTestTool"].tool
        sources = wf.step_nodes["CatTestTool"].sources
        process = translator.gen_process_from_cmdtool(tool, sources, scope=['CatTestTool'])
        expected_inputs = {'path inp'}
        actual_inputs = {inp.get_string() for inp in process.inputs}
        self.assertEqual(actual_inputs, expected_inputs)
    
    def test_static_inputs(self) -> None:
        # DO NOT need a process input for each static value in step sources.
        # non-files are fed data via params. 
        wf = AssemblyTestWF()
        refresh_workflow_inputs(wf)
        tool = wf.step_nodes["unicycler"].tool
        sources = wf.step_nodes["unicycler"].sources
        process = translator.gen_process_from_cmdtool(tool, sources, scope=['unicycler'])
        actual_inputs = {inp.name for inp in process.inputs}
        non_expected_ids = {
            'kmers',
            'scores',
            'startGeneCov',
            'startGeneId',
        }
        for input_id in non_expected_ids:
            self.assertNotIn(input_id, actual_inputs)
    
    def test_internal_inputs(self) -> None:
        # DO NOT need a process input for each non-exposed tool input. 
        # tool input will be autofilled or ignored in script.
        wf = AssemblyTestWF()
        refresh_workflow_inputs(wf)
        tool = wf.step_nodes["unicycler"].tool
        sources = wf.step_nodes["unicycler"].sources
        process = translator.gen_process_from_cmdtool(tool, sources, scope=['unicycler'])
        expected_ids = {
            'option1',
            'option2',
            'optionL',
        }
        actual_inputs = {inp.name for inp in process.inputs}
        non_expected_ids = {x.id() for x in tool.inputs() if x.id() not in expected_ids}
        for input_id in non_expected_ids:
            self.assertNotIn(input_id, actual_inputs)

    def test_arrays(self) -> None:
        # definition should be the same as singles. 
        # nextflow doesn't differentiate. 
        wf = ArrayStepInputsTestWF()
        refresh_workflow_inputs(wf)
        tool = wf.step_nodes["stp1"].tool
        sources = wf.step_nodes["stp1"].sources
        process = translator.gen_process_from_cmdtool(tool, sources, scope=['stp1'])
        actual_inputs = {inp.get_string() for inp in process.inputs}
        expected_inputs = {
            'path pos_basic',
            'path pos_basic2',
        }
        self.assertEqual(actual_inputs, expected_inputs)

    def test_secondaries(self) -> None:
        wf = SecondariesIOTestWF()
        refresh_workflow_inputs(wf)
        tool = wf.step_nodes["stp1"].tool
        sources = wf.step_nodes["stp1"].sources
        process = translator.gen_process_from_cmdtool(tool, sources, scope=['stp1'])
        actual_inputs = {inp.get_string() for inp in process.inputs}
        expected_inputs = {'tuple path(bam), path(bai)'}
        self.assertEqual(actual_inputs, expected_inputs)   
    
    def test_secondaries_array_inputs(self) -> None:
        wf = ArraySecondariesTestWF()
        refresh_workflow_inputs(wf)
        tool = wf.step_nodes["stp1"].tool
        sources = wf.step_nodes["stp1"].sources
        process = translator.gen_process_from_cmdtool(tool, sources, scope=['stp1'])
        actual_inputs = {inp.get_string() for inp in process.inputs}
        expected_inputs = {
            'path inp_bams',
            'path inp_bais'
        }
        self.assertEqual(actual_inputs, expected_inputs)   

    @unittest.skip('not implemented')
    def test_translate_commandtool(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_translate_codetool(self) -> None:
        raise NotImplementedError
    



class TestProcessOutputs(unittest.TestCase):
    """
    Need a process output for each tool output.
    """

    def setUp(self) -> None:
        settings.MINIMAL_PROCESS = True
        settings.MODE = 'workflow'
    
    def test_stdout(self):
        wf = BasicIOTestWF()
        refresh_workflow_inputs(wf)
        tool = wf.step_nodes["stp1"].tool
        sources = wf.step_nodes["stp1"].sources
        process = translator.gen_process_from_cmdtool(tool, sources, scope=['stp1'])
        actual_outputs = {inp.get_string() for inp in process.outputs}
        expected_outputs = {
            'stdout, emit: out',
        }
        self.assertEqual(actual_outputs, expected_outputs)  

    def test_wildcard_selector(self) -> None:
        wf = WildcardSelectorOutputTestWF()
        refresh_workflow_inputs(wf)
        # singles
        tool = wf.step_nodes["stp1"].tool
        sources = wf.step_nodes["stp1"].sources
        process = translator.gen_process_from_cmdtool(tool, sources, scope=['stp1'])
        actual_outputs = {out.get_string() for out in process.outputs}
        expected_outputs = {'path "myfile.txt", emit: out'}
        self.assertEqual(actual_outputs, expected_outputs)
        
        # arrays
        tool = wf.step_nodes["stp2"].tool
        sources = wf.step_nodes["stp2"].sources
        process = translator.gen_process_from_cmdtool(tool, sources, scope=['stp2'])
        actual_outputs = {out.get_string() for out in process.outputs}
        expected_outputs = {'path "*.txt", emit: out'}
        self.assertEqual(actual_outputs, expected_outputs)

    def test_input_selector(self) -> None:
        wf = InputSelectorOutputTestWF()
        refresh_workflow_inputs(wf)
        # singles
        tool = wf.step_nodes["stp1"].tool
        sources = wf.step_nodes["stp1"].sources
        process = translator.gen_process_from_cmdtool(tool, sources, scope=['stp1'])
        actual_outputs = {out.get_string() for out in process.outputs}
        expected_outputs = {'path inp, emit: out'}
        self.assertEqual(actual_outputs, expected_outputs)
        
        # arrays
        tool = wf.step_nodes["stp2"].tool
        sources = wf.step_nodes["stp2"].sources
        process = translator.gen_process_from_cmdtool(tool, sources, scope=['stp2'])
        actual_outputs = {out.get_string() for out in process.outputs}
        expected_outputs = {'path inp, emit: out'}
        self.assertEqual(actual_outputs, expected_outputs)

    def test_secondaries(self) -> None:
        wf = SecondariesIOTestWF()
        refresh_workflow_inputs(wf)
        tool = wf.step_nodes["stp1"].tool
        sources = wf.step_nodes["stp1"].sources
        process = translator.gen_process_from_cmdtool(tool, sources, scope=['stp1'])
        actual_outputs = {out.get_string() for out in process.outputs}
        expected_outputs = {'tuple path("*.bam"), path("*.bam.bai"), emit: out'}
        self.assertEqual(actual_outputs, expected_outputs)
    
    def test_secondaries_replaced(self) -> None:
        wf = SecondariesIOTestWF()
        refresh_workflow_inputs(wf)
        tool = wf.step_nodes["stp2"].tool
        sources = wf.step_nodes["stp2"].sources
        process = translator.gen_process_from_cmdtool(tool, sources, scope=['stp2'])
        actual_outputs = {out.get_string() for out in process.outputs}
        expected_outputs = {'tuple path("*.bam"), path("*.bai"), emit: out'}
        self.assertEqual(actual_outputs, expected_outputs)
    
    @unittest.skip('not implemented')
    def test_secondaries_array(self) -> None:
        # highly unlikely workflow would do this
        raise NotImplementedError

    @unittest.skip('not implemented')
    def test_translate_commandtool(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_translate_codetool(self) -> None:
        raise NotImplementedError




class TestProcessScript(unittest.TestCase):
    """
    Tests identifying tool inputs which should be process inputs.
    Need a process input for each tool input in step sources
    """

    def setUp(self) -> None:
        settings.MINIMAL_PROCESS = True
        settings.MODE = 'workflow'

    def test_filename_generated_tool(self):
        # TODO FIX ME
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

    @unittest.skip('not implemented')
    def test_translate_commandtool(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_translate_codetool(self) -> None:
        raise NotImplementedError




# distribute to the above TestProces___ classes
# class TestProcessCodeTool(unittest.TestCase):
#     def test_str_input(self):
#         p = translator.gen_process_from_codetool(SplitTextTestTool())
#         expected = f"""
# process TestSplitTextTool
# {{
#   input:
#     path PYTHON_CODE_FILE_PATH
#     val inp

#   output:
#     val "${{file("$workDir/janis_out_out").text.replace('[', '').replace(']', '').split(', ')}}" , emit: out

#   publishDir "$launchDir/TestSplitTextTool"
#   memory "${{params.runtime_memory ? params.runtime_memory + 'GB': ''}}"
#   cpus "${{params.runtime_cpu ? params.runtime_cpu : ''}}"
#   disk "${{params.runtime_disks ? params.runtime_disks + 'GB': ''}}"
#   time "${{params.runtime_seconds + 's'}}"

#   script:
#     \"\"\"
#     #!/usr/bin/env python
#     from TestSplitTextTool import code_block
#     import os
#     import json

#     result = code_block(inp="$inp")

#     work_dir = os.getenv("PYENV_DIR")
#     for key in result:
#         with open(os.path.join("$workDir", f"janis_out_{{key}}"), "w") as f:
#             f.write(json.dumps(result[key]))
#     \"\"\"
# }}


# """
#         self.assertEqual(expected, p.get_string())

#     def test_int_input(self):
#         p = translator.gen_process_from_codetool(SumTestTool())
#         expected = f"""
# process TestSumTool
# {{
#   input:
#     path PYTHON_CODE_FILE_PATH
#     val inp1
#     val inp2

#   output:
#     val "${{file("$workDir/janis_out_out").text.replace('[', '').replace(']', '').split(', ')}}" , emit: out

#   publishDir "$launchDir/TestSumTool"
#   memory "${{params.runtime_memory ? params.runtime_memory + 'GB': ''}}"
#   cpus "${{params.runtime_cpu ? params.runtime_cpu : ''}}"
#   disk "${{params.runtime_disks ? params.runtime_disks + 'GB': ''}}"
#   time "${{params.runtime_seconds + 's'}}"

#   script:
#     \"\"\"
#     #!/usr/bin/env python
#     from TestSumTool import code_block
#     import os
#     import json

#     result = code_block(inp1=$inp1, inp2=$inp2)

#     work_dir = os.getenv("PYENV_DIR")
#     for key in result:
#         with open(os.path.join("$workDir", f"janis_out_{{key}}"), "w") as f:
#             f.write(json.dumps(result[key]))
#     \"\"\"
# }}


# """
#         self.assertEqual(expected, p.get_string())

#     def test_array_input(self):
#         p = translator.gen_process_from_codetool(JoinArrayTestTool())
#         expected = f"""
# process TestJoinArrayTool
# {{
#   input:
#     path PYTHON_CODE_FILE_PATH
#     val inp

#   output:
#     val "${{file("$workDir/janis_out_out").text.replace('[', '').replace(']', '').split(', ')}}" , emit: out

#   publishDir "$launchDir/TestJoinArrayTool"
#   memory "${{params.runtime_memory ? params.runtime_memory + 'GB': ''}}"
#   cpus "${{params.runtime_cpu ? params.runtime_cpu : ''}}"
#   disk "${{params.runtime_disks ? params.runtime_disks + 'GB': ''}}"
#   time "${{params.runtime_seconds + 's'}}"

#   script:
#     \"\"\"
#     #!/usr/bin/env python
#     from TestJoinArrayTool import code_block
#     import os
#     import json

#     result = code_block(inp="$inp".split(" "))

#     work_dir = os.getenv("PYENV_DIR")
#     for key in result:
#         with open(os.path.join("$workDir", f"janis_out_{{key}}"), "w") as f:
#             f.write(json.dumps(result[key]))
#     \"\"\"
# }}


# """
#         self.assertEqual(expected, p.get_string())

#     def test_file_input(self):
#         p = translator.gen_process_from_codetool(FileInputTestTool())
#         expected = f"""
# process TestFileInput
# {{
#   input:
#     path PYTHON_CODE_FILE_PATH
#     path inp

#   output:
#     val "${{file("$workDir/janis_out_out").text.replace('[', '').replace(']', '').split(', ')}}" , emit: out

#   publishDir "$launchDir/TestFileInput"
#   memory "${{params.runtime_memory ? params.runtime_memory + 'GB': ''}}"
#   cpus "${{params.runtime_cpu ? params.runtime_cpu : ''}}"
#   disk "${{params.runtime_disks ? params.runtime_disks + 'GB': ''}}"
#   time "${{params.runtime_seconds + 's'}}"

#   script:
#     \"\"\"
#     #!/usr/bin/env python
#     from TestFileInput import code_block
#     import os
#     import json

#     result = code_block(inp="$inp")

#     work_dir = os.getenv("PYENV_DIR")
#     for key in result:
#         with open(os.path.join("$workDir", f"janis_out_{{key}}"), "w") as f:
#             f.write(json.dumps(result[key]))
#     \"\"\"
# }}


# """
#         self.assertEqual(expected, p.get_string())

#     def test_file_with_secondary_input(self):
#         p = translator.gen_process_from_codetool(
#             SecondaryInputTestTool()
#         )
#         expected = f"""
# process TestFileWithSecondaryInput
# {{
#   input:
#     path PYTHON_CODE_FILE_PATH
#     path inp

#   output:
#     val "${{file("$workDir/janis_out_out").text.replace('[', '').replace(']', '').split(', ')}}" , emit: out

#   publishDir "$launchDir/TestFileWithSecondaryInput"
#   memory "${{params.runtime_memory ? params.runtime_memory + 'GB': ''}}"
#   cpus "${{params.runtime_cpu ? params.runtime_cpu : ''}}"
#   disk "${{params.runtime_disks ? params.runtime_disks + 'GB': ''}}"
#   time "${{params.runtime_seconds + 's'}}"

#   script:
#     \"\"\"
#     #!/usr/bin/env python
#     from TestFileWithSecondaryInput import code_block
#     import os
#     import json

#     result = code_block(inp="$inp".split(" ")[0])

#     work_dir = os.getenv("PYENV_DIR")
#     for key in result:
#         with open(os.path.join("$workDir", f"janis_out_{{key}}"), "w") as f:
#             f.write(json.dumps(result[key]))
#     \"\"\"
# }}


# """
#         self.assertEqual(expected, p.get_string())




class TestPlumbingBasic(unittest.TestCase):
    """
    This test group checks janis 'TInput' step inputs driven by workflow inputs
    or static values. 

    Ensures they are being handled correctly when parsed to nextflow.
    """

    def setUp(self) -> None:
        settings.MINIMAL_PROCESS = True
        settings.MODE = 'workflow'

    # workflow input step inputs
    def test_workflow_inputs(self):
        wf = StepInputsWFInputTestWF()
        refresh_workflow_inputs(wf)
        tool = wf.step_nodes["stp1"].tool
        sources = wf.step_nodes["stp1"].sources
        scatter = wf.step_nodes["stp1"].scatter
        expected = {
            "pos_basic": "ch_in_file",
            "pos_basic2": "ch_in_file_opt",
        }
        actual = translator.gen_step_inval_dict(tool, sources, scatter)
        self.assertEqual(expected, actual)
        
    # static step inputs
    def test_static_inputs(self):
        wf = StepInputsTestWF()
        refresh_workflow_inputs(wf)
        tool = wf.step_nodes["stp2"].tool
        sources = wf.step_nodes["stp2"].sources
        scatter = wf.step_nodes["stp2"].scatter
        not_expected = {
            'pos_default',
            'pos_default',
            'pos_optional',
            'flag_true',
            'flag_false',
            'opt_basic',
            'opt_default',
            'opt_optional',
        }
        actual = translator.gen_step_inval_dict(tool, sources, scatter)
        for tinput_name in not_expected:
            self.assertNotIn(tinput_name, actual)

    # connections
    def test_connections_files(self) -> None:
        wf = StepConnectionsTestWF()
        refresh_workflow_inputs(wf)
        tool = wf.step_nodes["stp2"].tool
        sources = wf.step_nodes["stp2"].sources
        scatter = wf.step_nodes["stp2"].scatter
        expected = {"inp": "STP1.out.out"}
        actual = translator.gen_step_inval_dict(tool, sources, scatter)
        self.assertEqual(expected, actual)

    @unittest.skip('not implemented')
    def test_connections_nonfiles(self) -> None:
        raise NotImplementedError



class TestPlumbingBasicArrays(unittest.TestCase):
    """
    This test group checks janis 'TInput' step inputs driven by workflow inputs
    or static values as above, but for Array situations. 

    Ensures they are being handled correctly when parsed to nextflow.
    """

    def setUp(self) -> None:
        settings.MINIMAL_PROCESS = True
        settings.MODE = 'workflow'

    def test_array_connections(self) -> None:
        wf = ArrayStepConnectionsTestWF()
        refresh_workflow_inputs(wf)
        tool = wf.step_nodes["stp2"].tool
        sources = wf.step_nodes["stp2"].sources
        scatter = wf.step_nodes["stp2"].scatter
        expected = {"inp": "STP1.out.out"}
        actual = translator.gen_step_inval_dict(tool, sources, scatter)
        self.assertEqual(expected, actual)
    
    def test_workflow_inputs_array(self) -> None:
        wf = ArrayStepInputsTestWF()
        refresh_workflow_inputs(wf)
        tool = wf.step_nodes["stp1"].tool
        sources = wf.step_nodes["stp1"].sources
        scatter = wf.step_nodes["stp1"].scatter
        expected = {
            "pos_basic": "ch_in_file_array",
            "pos_basic2": "ch_in_file_array_opt",
        }
        actual = translator.gen_step_inval_dict(tool, sources, scatter)
        self.assertEqual(expected, actual)

    def test_static_step_inputs_array(self):
        wf = ArrayStepInputsTestWF()
        refresh_workflow_inputs(wf)
        tool = wf.step_nodes["stp2"].tool
        sources = wf.step_nodes["stp2"].sources
        scatter = wf.step_nodes["stp2"].scatter
        not_expected = {
            'pos_default',
            'pos_default',
            'pos_optional',
            'flag_true',
            'flag_false',
            'opt_basic',
            'opt_default',
            'opt_optional',
        }
        actual = translator.gen_step_inval_dict(tool, sources, scatter)
        for tinput_name in not_expected:
            self.assertNotIn(tinput_name, actual)


class TestPlumbingScatter(unittest.TestCase):
    """
    This test group checks that janis scatter declarations are being handled 
    correctly to produce the desired nextflow workflow. 
    """

    def setUp(self) -> None:
        settings.MINIMAL_PROCESS = True
        settings.MODE = 'workflow'
    
    def test_scatter_basic(self) -> None:
        wf = BasicScatterTestWF()
        refresh_workflow_inputs(wf)
        tool = wf.step_nodes["stp1"].tool
        sources = wf.step_nodes["stp1"].sources
        scatter = wf.step_nodes["stp1"].scatter
        expected = {"inp": "ch_in_file_array.flatten()"}
        actual = translator.gen_step_inval_dict(tool, sources, scatter)
        self.assertEqual(expected, actual)
    
    def test_scatter_dot(self) -> None:
        wf = ScatterDotTestWF()
        refresh_workflow_inputs(wf)
        tool = wf.step_nodes["stp1"].tool
        sources = wf.step_nodes["stp1"].sources
        scatter = wf.step_nodes["stp1"].scatter
        expected = {
            "pos_basic": "ch_in_file_array.flatten()",
            "opt_basic": "ch_in_str_array.flatten()",
        }
        actual = translator.gen_step_inval_dict(tool, sources, scatter)
        self.assertEqual(expected, actual)
    
    def test_scatter_cross(self) -> None:
        wf = ScatterCrossTestWF()
        refresh_workflow_inputs(wf)
        tool = wf.step_nodes["stp1"].tool
        sources = wf.step_nodes["stp1"].sources
        scatter = wf.step_nodes["stp1"].scatter

        # cartesian cross channel manipulation in workflow
        operation = nfgen.channels.gen_scatter_cross_operation(sources, scatter)
        actual_op = operation.get_string()
        expected_op = """\
ch_in_file_array.flatten()
.combine(ch_in_str_array).flatten()
.multiMap { it ->
    in_file_array: it[0]
    in_str_array: it[1]
}
.set { ch_cartesian_cross }"""
        self.assertEqual(expected_op, actual_op)

        # step input values
        expected = {
            "pos_basic": "ch_cartesian_cross.in_file_array",
            "opt_basic": "ch_cartesian_cross.in_str_array",
        }
        actual = translator.gen_step_inval_dict(tool, sources, scatter)
        self.assertEqual(expected, actual)
    
    def test_scatter_connection(self) -> None:
        wf = ChainedScatterTestWF()
        refresh_workflow_inputs(wf)
        # single -> single
        tool = wf.step_nodes["stp2"].tool
        sources = wf.step_nodes["stp2"].sources
        scatter = wf.step_nodes["stp2"].scatter
        expected = {"inp": "STP1.out.out"}
        actual = translator.gen_step_inval_dict(tool, sources, scatter)
        self.assertEqual(expected, actual)
        # array -> single
        tool = wf.step_nodes["stp4"].tool
        sources = wf.step_nodes["stp4"].sources
        scatter = wf.step_nodes["stp4"].scatter
        expected = {"inp": "STP3.out.out.flatten()"}
        actual = translator.gen_step_inval_dict(tool, sources, scatter)
        self.assertEqual(expected, actual)




class TestPlumbingSecondaries(unittest.TestCase):
    # TODO test more features / edge cases?
    """
    This test group checks that janis File types with secondaries are being handled 
    correctly to produce the desired nextflow workflow. 
    
    SecondaryFiles are especially tricky when mapping to janis. 
    Requires the remapping of tool inputs in some situations.
    
    eg Array(BamBai) 
        - Arrays where each item is a file with secondaries. 
        - These gets remapped from a single tool input, to an array input
          for each separate file in the datatype. 
        - Works like transposing - instead of inp1=[[Bam, Bai], [Bam, Bai]],
          results in inp1=[Bam, Bam], inp2=[Bai, Bai]

    """

    def setUp(self) -> None:
        settings.MINIMAL_PROCESS = True
        settings.MODE = 'workflow'
    
    def test_secondaries_workflow_inputs(self) -> None:
        wf = SecondariesIOTestWF()
        refresh_workflow_inputs(wf)
        # params created correctly?
        actual_params = nfgen.params.serialize()
        expected_params = {
            'in_alignments_bam': 'null',
            'in_alignments_bai': 'null',
        }
        self.assertEquals(actual_params, expected_params)
        
        # channel created correctly?
        reference = 'inAlignments'
        expected_channel = 'ch_in_alignments = Channel.fromPath( params.in_alignments_bam, params.in_alignments_bai ).collect()'
        channel = nfgen.channels.getall(reference)[0]
        actual_channel = channel.declaration
        self.assertEquals(actual_channel, expected_channel)
        
        # step inputs created correctly?
        tool = wf.step_nodes["stp1"].tool
        sources = wf.step_nodes["stp1"].sources
        scatter = wf.step_nodes["stp1"].scatter
        expected_inputs = {"inp": "ch_in_alignments"}
        actual_inputs = translator.gen_step_inval_dict(tool, sources, scatter)
        self.assertEquals(actual_inputs, expected_inputs)
    
    def test_secondaries_connections(self) -> None:
        wf = SecondariesConnectionsTestWF()
        refresh_workflow_inputs(wf)
        # step inputs created correctly?
        tool = wf.step_nodes["stp2"].tool
        sources = wf.step_nodes["stp2"].sources
        scatter = wf.step_nodes["stp2"].scatter
        expected_inputs = {"inp": "STP1.out.out"}
        actual_inputs = translator.gen_step_inval_dict(tool, sources, scatter)
        self.assertEquals(actual_inputs, expected_inputs)




class TestPlumbingCombinations(unittest.TestCase):
    """
    Tests plumbing for most complex cases.
    Includes combinations of scatter, arrays, secondary files. 
    """
    def setUp(self) -> None:
        settings.MINIMAL_PROCESS = True
        settings.MODE = 'workflow'

    # scatter + secondaries
    @unittest.skip('not implemented')
    def test_scatter_secondaries_workflow_inputs(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_scatter_secondaries_dot(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_scatter_secondaries_cross(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_scatter_secondaries_connections(self) -> None:
        raise NotImplementedError

    # secondaries + arrays
    @unittest.skip('not implemented')
    def test_secondaries_workflow_inputs_array(self) -> None:
        # wf = ArraySecondariesTestWF()
        # actual_params = refresh_workflow_inputs(wf)
        # expected_params = {
        #     'in_alignments_bams': '[]',
        #     'in_alignments_bais': '[]',
        # }
        # self.assertEquals(actual_params, expected_params)
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_secondaries_connections_array(self) -> None:
        raise NotImplementedError
    
    # scatter + arrays
    @unittest.skip('not implemented')
    def test_scatter_basic_array(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_scatter_dot_array(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_scatter_cross_array(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_scatter_connection_array(self) -> None:
        raise NotImplementedError



class TestWorkflowOutputs(unittest.TestCase):
    """
    Tests workflow outputs being created correctly
    """
    def setUp(self) -> None:
        settings.MINIMAL_PROCESS = True
        settings.MODE = 'workflow'
        
    @unittest.skip('not implemented')
    def test_files(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_files_arrays(self) -> None:
        raise NotImplementedError

    @unittest.skip('not implemented')
    def test_secondaries(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_secondaries_array(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_scatter(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_scatter_secondaries(self) -> None:
        raise NotImplementedError




class TestNextflowConfig(unittest.TestCase):
    
    def setUp(self) -> None:
        settings.MINIMAL_PROCESS = True
        settings.MODE = 'workflow'
    
    def test_file_workflow_inputs(self):
        wf = BasicIOTestWF()
        refresh_workflow_inputs(wf)
        config = translator.stringify_translated_inputs({})
        expected_values = {
            'in_file': 'null',
            'in_str': 'null',
            'in_int': 'null',
            'stp3_inp': '50',
        }
        for name, val in expected_values.items():
            pattern = f'{name}.*?{val}'
            matches = re.findall(pattern, config)
            self.assertGreater(len(matches), 0)
    
    def test_array_workflow_inputs(self):
        wf = ArrayIOExtrasTestWF()
        refresh_workflow_inputs(wf)
        config = translator.stringify_translated_inputs({})
        # I apologise for disgusting formatting below. 
        # This can probably be handled better.
        expected_values = {
    'in_str_array': fr"""\[
        // list strings here
    \]""",
    'in_file_array': fr"""\[
        // list files here
    \]""",
            'in_int_array': r'\[\]',
            'in_float_array': r'\[\]',
        }
        print(config)
        for name, val in expected_values.items():
            pattern = f'{name}.*?{val}'
            matches = re.findall(pattern, config)
            self.assertGreater(len(matches), 0)
    
    def test_file_secondaries_workflow_inputs(self):
        wf = SecondariesIOTestWF()
        refresh_workflow_inputs(wf)
        config = translator.stringify_translated_inputs({})
        expected_values = {
            'in_alignments_bam': 'null',
            'in_alignments_bai': 'null',
        }
        print(config)
        for name, val in expected_values.items():
            pattern = f'{name}.*?{val}'
            matches = re.findall(pattern, config)
            self.assertGreater(len(matches), 0)
    
    def test_file_secondaries_array_workflow_inputs(self):
        wf = ArraySecondariesTestWF()
        refresh_workflow_inputs(wf)
        config = translator.stringify_translated_inputs({})
        expected_values = {
    'in_alignments_bams': fr"""\[
        // list files here
    \]""",
    'in_alignments_bais': fr"""\[
        // list files here
    \]""",
        }
        print(config)
        for name, val in expected_values.items():
            pattern = f'{name}.*?{val}'
            matches = re.findall(pattern, config)
            self.assertGreater(len(matches), 0)
    
    def test_nonfile_workflow_inputs(self):
        # string, int, bool
        wf = StepInputsTestWF()
        refresh_workflow_inputs(wf)
        config = translator.stringify_translated_inputs({})
        expected_values = {
            'in_file': 'null',
            'in_str': 'null',
            'in_int': 'null',
            'in_bool': 'null',
            'stp2_pos_default': '100',
            'stp2_pos_optional': "'static'",
            'stp2_flag_true': 'false',
            'stp2_flag_false': 'true',
            'stp2_opt_basic': "'static'",
            'stp2_opt_default': '100',
            'stp2_opt_optional': "''",
            'stp3_opt_basic': "'static'",
            'stp3_opt_default': '100',
            'stp3_opt_optional': "''",
        }
        print(config)
        for name, val in expected_values.items():
            pattern = f'{name}.*?{val}'
            matches = re.findall(pattern, config)
            self.assertGreater(len(matches), 0)
    
    def test_nonfile_array_workflow_inputs(self):
        # string, int, bool
        wf = ArrayStepInputsTestWF()
        refresh_workflow_inputs(wf)
        config = translator.stringify_translated_inputs({})
        expected_values = {
    'in_file_array': fr"""\[
        // list files here
    \]""",
    'in_str_array': fr"""\[
        // list strings here
    \]""",
    'in_int_array': fr'\[\]',
    'in_bool_array': fr'\[\]',
    'stp2_pos_default': fr'\[4, 5, 6\]',
    'stp2_pos_optional': fr"""\[
        'hi',
        'there',
        'friend',
    \]""",
    'stp2_flag_true': fr'\[true\]',
    'stp2_flag_false': fr'\[true\]',
    'stp2_opt_basic': fr"""\[
        'hi',
        'there',
        'friend',
    \]""",
    'stp2_opt_default': fr'\[4, 5, 6\]',
    'stp2_opt_optional': fr"""\[
        'hi',
        'there',
        'friend',
    \]""",
        }
        print(config)
        for name, val in expected_values.items():
            pattern = f'{name}.*?{val}'
            matches = re.findall(pattern, config)
            self.assertGreater(len(matches), 0)

    def test_workflow_config(self) -> None:
        wf = AssemblyTestWF()
        refresh_workflow_inputs(wf)
        params = translator.build_inputs_file(wf)
        config = translator.stringify_translated_inputs(params)
        # basic structure
        self.assertIn('docker.enabled = true', config)
        self.assertIn('params {\n\n', config)
        self.assertIn('\n\n}', config)
        # expected params are present & correct
        expected_values = {
            'test_input': 'null',
            'fastqc1_adapters': 'null',
            'fastqc1_contaminants': 'null',
            'fastqc1_limits': 'null',
            'fastqc2_adapters': 'null',
            'fastqc2_contaminants': 'null',
            'fastqc2_limits': 'null',
            'in_forward_reads': 'null',
            'in_long_reads': 'null',
            'in_reverse_reads': 'null',
            'unicycler_kmers': "''",
            'unicycler_scores': "''",
            'unicycler_start_gene_cov': '95.0',
            'unicycler_start_gene_id': '90.0',
        }
        for name, val in expected_values.items():
            pattern = f'{name}.*?{val}'
            matches = re.findall(pattern, config)
            self.assertGreater(len(matches), 0)

    @unittest.skip('not implemented')
    def test_tool_config(self) -> None:
        settings.MINIMAL_PROCESS = True
        settings.MODE = 'tool'
        raise NotImplementedError




class TestWorkflow(unittest.TestCase):
    def setUp(self) -> None:
        settings.MINIMAL_PROCESS = True
        settings.MODE = 'workflow'




class TestSubworkflows(unittest.TestCase):
    # TODO test more features / edge cases?

    def setUp(self) -> None:
        settings.MINIMAL_PROCESS = True
        settings.MODE = 'workflow'




class TestStepFeatures(unittest.TestCase):

    def setUp(self) -> None:
        settings.MINIMAL_PROCESS = True
        settings.MODE = 'workflow'

    def test_first_selector(self):
        workflow = ConditionStepTestWF()

        step_id = "print"
        tool = workflow.step_nodes[step_id].tool
        sources = workflow.step_nodes[step_id].sources
        scatter = workflow.step_nodes[step_id].scatter
        inputs = translator.gen_step_inval_dict(tool, sources, scatter)
        expected = {"inp": "[$params.mystring, $get_string.out.out].first()"}

        self.assertEqual(expected, inputs)

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




class TestUnwrap(unittest.TestCase):
    any_tool = BasicTestTool()

    def test_string_formatter(self):
        b = StringFormatter("no format")
        res = nfgen.translate_string_formatter(b, self.any_tool, {})
        self.assertEqual("no format", res)

    def test_string_formatter_one_string_param(self):
        b = StringFormatter("there's {one} arg", one="a string")
        res = nfgen.translate_string_formatter(b, self.any_tool, {})
        self.assertEqual("there's ${'a string'} arg", res)

    def test_string_formatter_one_input_selector_param(self):
        b = StringFormatter("an input {arg}", arg=InputSelector("testtool"))
        res = nfgen.translate_string_formatter(
            b, self.any_tool, input_in_selectors={}, inputs_dict=self.any_tool.inputs_map(),
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
            b, tool, input_in_selectors={}, inputs_dict=tool.inputs_map()
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
            b, tool, input_in_selectors={}, inputs_dict=tool.inputs_map()
        )
        self.assertEqual("${user}\\t${static}", res)

        res2 = nfgen.translate_string_formatter(
            b, tool, input_in_selectors={}, inputs_dict=tool.inputs_map(),
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
            b, tool, input_in_selectors={}, inputs_dict=tool.inputs_map()
        )
        self.assertEqual("${testtool}:${arrayInp.join(';')}", res)




"""
FROM TestNextflowConfig

    def test_input_in_input_value_nooptional_default(self):
        wf = WorkflowBuilder("test_nf_input_in_input_value_nooptional_default")
        wf.input("inpId", String(), default="1")
        ai = {"inpId": "2"}
        self.assertDictEqual(
            {"inpId": "2"}, translator.build_inputs_file(wf, additional_inputs=ai)
        )

    def test_input_in_input_value_optional_nodefault(self):
        wf = WorkflowBuilder("test_nf_input_in_input_value_optional_nodefault")
        wf.input("inpId", String(optional=True))
        ai = {"inpId": "2"}
        self.assertDictEqual(
            {"inpId": "2"}, translator.build_inputs_file(wf, additional_inputs=ai)
        )

    def test_input_in_input_value_optional_default(self):
        wf = WorkflowBuilder("test_nf_input_in_input_value_optional_default")
        wf.input("inpId", String(optional=True), default="1")
        ai = {"inpId": "2"}
        self.assertDictEqual(
            {"inpId": "2"}, translator.build_inputs_file(wf, additional_inputs=ai)
        )

    def test_input_in_input_novalue_nooptional_nodefault(self):
        wf = WorkflowBuilder("test_nf_input_in_input_novalue_nooptional_nodefault")
        wf.input("inpId", String())
        # included because no value, no default, and not optional
        self.assertDictEqual({"inpId": ""}, translator.build_inputs_file(wf))

    def test_input_in_input_novalue_nooptional_default(self):
        wf = WorkflowBuilder("test_nf_input_in_input_novalue_nooptional_default")
        wf.input("inpId", String(), default="1")
        self.assertDictEqual({"inpId": "1"}, translator.build_inputs_file(wf))

    def test_input_in_input_novalue_optional_nodefault(self):
        wf = WorkflowBuilder("test_nf_input_in_input_novalue_optional_nodefault")
        wf.input("inpId", String(optional=True))
        self.assertDictEqual({"inpId": ""}, translator.build_inputs_file(wf))

    def test_input_in_input_novalue_optional_default(self):
        wf = WorkflowBuilder("test_nf_input_in_input_novalue_optional_default")
        wf.input("inpId", String(optional=True), default="1")
        self.assertDictEqual({"inpId": "1"}, translator.build_inputs_file(wf))

    def test_bool_input_type_value(self):
        wf = WorkflowBuilder("test_bool_input")
        wf.input("inpBool", Boolean())
        ai = {"inpBool": True}
        self.assertDictEqual(
            ai, translator.build_inputs_file(wf, additional_inputs=ai)
        )

    def test_bool_input_type_novalue_default(self):
        wf = WorkflowBuilder("test_bool_input")
        wf.input("inpBool", Boolean(), default=False)
        self.assertDictEqual({"inpBool": False}, translator.build_inputs_file(wf))

    def test_bool_input_type_value_default(self):
        wf = WorkflowBuilder("test_bool_input")
        wf.input("inpBool", Boolean(), default=False)
        ai = {"inpBool": True}
        self.assertDictEqual(
            ai, translator.build_inputs_file(wf, additional_inputs=ai)
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
            translator.build_inputs_file(wf, additional_inputs=ai),
        )

    def test_file_type_no_secondary_value(self):
        wf = WorkflowBuilder("test_file_input")
        wf.input("inp", DataTypeNoSecondary())
        ai = {"inp": "/my/path/filename.doc"}
        self.assertDictEqual(
            {"inp": "/my/path/filename.doc"},
            translator.build_inputs_file(wf, additional_inputs=ai),
        )

    def test_file_type_with_secondary_no_value(self):
        wf = WorkflowBuilder("test_file_input")
        wf.input("inp", DataTypeWithSecondary())
        self.assertDictEqual(
            {"inp": f"/{settings.NO_FILE_PATH_PREFIX}1"},
            translator.build_inputs_file(wf),
        )

    def test_file_type_no_secondary_no_value(self):
        wf = WorkflowBuilder("test_file_input")
        wf.input("inp", DataTypeNoSecondary())
        self.assertDictEqual(
            {"inp": f"/{settings.NO_FILE_PATH_PREFIX}1"},
            translator.build_inputs_file(wf),
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
            translator.build_inputs_file(wf),
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
            translator.build_inputs_file(wf),



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




class TestPrepareInputVars(unittest.TestCase):
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
"""