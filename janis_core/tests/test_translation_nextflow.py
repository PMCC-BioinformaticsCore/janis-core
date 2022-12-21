import unittest
import regex as re

from janis_core.tests.testtools import (
    InputQualityTestTool,
    BasicTestTool,
)

from janis_core.tests.testworkflows import (

    # basics
    BasicIOTestWF,
    WildcardSelectorOutputTestWF,
    InputSelectorTestWF,
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
    IndexOperatorTestWF,
    StringFormatterTestWF,
    FilenameGeneratedTestWF,

    # codetools
    InputsPythonToolTestWF,
    OutputsPythonToolTestWF,

    # specific workflows
    AssemblyTestWF,
    SubworkflowTestWF,
    FilenameTestWF,
    OutputCollectionTestWF,
    UnwrapTestWF,
)

from janis_core import (
    Workflow,
    InputSelector,
    StringFormatter,
    JoinOperator,
)

from janis_core.translations import NextflowTranslator as translator
from janis_core.translations import nfgen
from janis_core.translations.nfgen import settings
from janis_core import (
    String, 
    Int, 
    Float, 
    File, 
    Boolean, 
)
from janis_core.translations.nfgen.nfgen_utils import to_groovy


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
    scope = [settings.NF_MAIN_NAME]
    nfgen.params.clear()
    nfgen.channels.clear()
    nfgen.register_params_channels(wf, scope=scope)



class TestToGroovyStr(unittest.TestCase):

    def test_string(self) -> None:
        inp = 'Hello'
        expected = "'Hello'"
        actual = to_groovy(inp, String())
        self.assertEquals(expected, actual)
    
    def test_numeric(self) -> None:
        # int
        inp = 5
        expected = '5'
        actual = to_groovy(inp, Int())
        self.assertEquals(expected, actual)
        # float
        inp = 5.2
        expected = '5.2'
        actual = to_groovy(inp, Float())
        self.assertEquals(expected, actual)
    
    def test_bool(self) -> None:
        # true
        inp = True
        expected = 'true'
        actual = to_groovy(inp, Boolean())
        self.assertEquals(expected, actual)
        # false
        inp = False
        expected = 'false'
        actual = to_groovy(inp, Boolean())
        self.assertEquals(expected, actual)
    
    def test_none(self) -> None:
        inp = None
        expected = 'null'
        actual = to_groovy(inp, String())
        self.assertEquals(expected, actual)
    
    def test_empty_array(self) -> None:
        inp = []
        expected = "[]"
        actual = to_groovy(inp, String())
        self.assertEquals(expected, actual)
    
    def test_string_array(self) -> None:
        inp = ['Hello']
        expected = "['Hello']"
        actual = to_groovy(inp, String())
        self.assertEquals(expected, actual)
    
    def test_numeric_array(self) -> None:
        # int
        inp = [1, 2, 3]
        expected = '[1, 2, 3]'
        actual = to_groovy(inp, Int())
        self.assertEquals(expected, actual)
        # float
        inp = [1.1, 2.2, 3.3]
        expected = '[1.1, 2.2, 3.3]'
        actual = to_groovy(inp, Float())
        self.assertEquals(expected, actual)
    
    def test_bool_array(self) -> None:
        inp = [True, False]
        expected = '[true, false]'
        actual = to_groovy(inp)
        self.assertEquals(expected, actual)
    
    def test_none_array(self) -> None:
        inp = [None, None]
        expected = '[null, null]'
        actual = to_groovy(inp)
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
    def test_translate_pythontool(self) -> None:
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
        expected_channel_names = {
            'ch_fastqc1_adapters',
            'ch_fastqc1_contaminants',
            'ch_fastqc1_limits',
            'ch_fastqc2_adapters',
            'ch_fastqc2_contaminants',
            'ch_fastqc2_limits',
        }
        all_channels = nfgen.channels.getall()
        optional_channels = [ch for ch in all_channels if ch.name in expected_channel_names]

        # check each expected channel exists
        self.assertEqual(len(optional_channels), len(expected_channel_names))
        
        # for each optional channel, check it has correct format
        for channel in optional_channels:
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
        inp = wf.input_nodes['inAlignments']
        alignments_ch = nfgen.channels.get(inp.uuid)
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

    def test_filename_types(self) -> None:
        wf = FilenameTestWF()
        refresh_workflow_inputs(wf)
        channels_ids = {c.name for c in nfgen.channels.getall()}
        expected_ids = {
            'ch_in_file',
            'ch_in_str',
        }
        self.assertEqual(channels_ids, expected_ids)

    @unittest.skip('not implemented')
    def test_channel_methods(self) -> None:
        raise NotImplementedError

    @unittest.skip('not implemented')
    def test_translate_commandtool(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_translate_pythontool(self) -> None:
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
        directives = nfgen.ordering.order_nf_directives(process.directives)
        actual_order = [type(x).__name__ for x in directives]
        expected_order = [
            'ContainerDirective',
            'PublishDirDirective',
            'DebugDirective',
            'CpusDirective',
            'DiskDirective',
            'MemoryDirective',
            'TimeDirective',
        ]
        for actual, expected in zip(actual_order, expected_order):
            self.assertEqual(actual, expected)

    def test_directives(self) -> None:
        wf = DirectivesTestWF()
        refresh_workflow_inputs(wf)
        step_id = 'stp1'
        tool = wf.step_nodes[step_id].tool
        sources = wf.step_nodes[step_id].sources
        process = translator.gen_process_from_cmdtool(tool, sources, scope=[settings.NF_MAIN_NAME, step_id])
        process = translator.handle_container(tool, process)
        actual_directives = {d.get_string() for d in process.directives}
        expected_directives = {
            'container "quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0"',
            'publishDir "${params.outdir}/stp1"',
            'debug true',
            'disk "${params.stp1.disk}"',
            'memory "${params.stp1.memory}"',
            'time "${params.stp1.time}"'
        }
        for direc in expected_directives:
            self.assertIn(direc, actual_directives)
    
    @unittest.skip('not implemented')
    def test_translate_commandtool(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_translate_pythontool(self) -> None:
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
    
    def test_pythontool(self) -> None:
        wf = InputsPythonToolTestWF()
        refresh_workflow_inputs(wf)
        
        # File, String, Int input types
        tool = wf.step_nodes["stp0"].tool
        sources = wf.step_nodes["stp0"].sources
        process = translator.gen_process_from_codetool(tool, sources, scope=['stp0'])
        actual_inputs = {inp.get_string() for inp in process.inputs}
        expected_inputs = {
            'path PYTHON_CODE_FILE_PATH',
            'path inp1',
        }
        self.assertEqual(actual_inputs, expected_inputs)
        
        # Array(String) input type
        tool = wf.step_nodes["stp1"].tool
        sources = wf.step_nodes["stp1"].sources
        process = translator.gen_process_from_codetool(tool, sources, scope=['stp1'])
        actual_inputs = {inp.get_string() for inp in process.inputs}
        expected_inputs = {
            'path PYTHON_CODE_FILE_PATH',
        }
        self.assertEqual(actual_inputs, expected_inputs)
        
        # File (secondaries) input type
        tool = wf.step_nodes["stp2"].tool
        sources = wf.step_nodes["stp2"].sources
        process = translator.gen_process_from_codetool(tool, sources, scope=['stp2'])
        actual_inputs = {inp.get_string() for inp in process.inputs}
        expected_inputs = {
            'path PYTHON_CODE_FILE_PATH',
            'tuple path(bam), path(bai)',
        }
        self.assertEqual(actual_inputs, expected_inputs)

    def test_filename_types(self) -> None:
        wf = FilenameTestWF()
        refresh_workflow_inputs(wf)

        step_id = 'stp1'
        tool = wf.step_nodes[step_id].tool
        sources = wf.step_nodes[step_id].sources
        process = translator.gen_process_from_cmdtool(tool, sources, scope=[settings.NF_MAIN_NAME, step_id])
        actual_inputs = {inp.get_string() for inp in process.inputs}
        expected_inputs = {
            'path inp1',
            'val inp2',
        }
        self.assertEqual(actual_inputs, expected_inputs)

        step_id = 'stp2'
        tool = wf.step_nodes[step_id].tool
        sources = wf.step_nodes[step_id].sources
        process = translator.gen_process_from_cmdtool(tool, sources, scope=[settings.NF_MAIN_NAME, step_id])
        actual_inputs = {inp.get_string() for inp in process.inputs}
        expected_inputs = {
            'path inp1',
        }
        self.assertEqual(actual_inputs, expected_inputs)

    @unittest.skip('not implemented')
    def test_translate_commandtool(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_translate_pythontool(self) -> None:
        raise NotImplementedError
    



class TestProcessOutputs(unittest.TestCase):
    """
    Need a process output for each tool output.
    """

    def setUp(self) -> None:
        settings.MINIMAL_PROCESS = True
        settings.MODE = 'workflow'
        self.wf = OutputCollectionTestWF()
        refresh_workflow_inputs(self.wf)

    def test_stdout(self):
        wf = BasicIOTestWF()
        refresh_workflow_inputs(wf)
        step_id = "stp1"
        tool = wf.step_nodes[step_id].tool
        sources = wf.step_nodes[step_id].sources
        scope = [settings.NF_MAIN_NAME, step_id]
        process = translator.gen_process_from_cmdtool(tool, sources, scope)
        actual_outputs = {out.get_string() for out in process.outputs}
        expected_outputs = {'stdout, emit: out'}
        self.assertEqual(actual_outputs, expected_outputs)  

    def test_wildcard(self) -> None:
        step_id = "stp1"
        tool = self.wf.step_nodes[step_id].tool
        sources = self.wf.step_nodes[step_id].sources
        scope = [settings.NF_MAIN_NAME, step_id]
        process = translator.gen_process_from_cmdtool(tool, sources, scope)
        actual_outputs = {out.get_string() for out in process.outputs}
        expected_outputs = {'path "myfile.txt", emit: out'}
        self.assertEqual(actual_outputs, expected_outputs)
    
    def test_wildcard_array(self) -> None:
        wf = WildcardSelectorOutputTestWF()
        refresh_workflow_inputs(wf)
        tool = wf.step_nodes["stp2"].tool
        sources = wf.step_nodes["stp2"].sources
        process = translator.gen_process_from_cmdtool(tool, sources, scope=['stp2'])
        actual_outputs = {out.get_string() for out in process.outputs}
        expected_outputs = {'path "*.txt", emit: out'}
        self.assertEqual(actual_outputs, expected_outputs)
    
    def test_input_selector(self) -> None:
        wf = InputSelectorTestWF()
        refresh_workflow_inputs(wf)
        tool = wf.step_nodes["stp1"].tool
        sources = wf.step_nodes["stp1"].sources
        process = translator.gen_process_from_cmdtool(tool, sources, scope=['stp1'])
        actual_outputs = {out.get_string() for out in process.outputs}
        expected_outputs = {'path inp, emit: out'}
        self.assertEqual(actual_outputs, expected_outputs)

    def test_input_selector_param(self) -> None:
        step_id = "stp4"
        tool = self.wf.step_nodes[step_id].tool
        sources = self.wf.step_nodes[step_id].sources
        scope = [settings.NF_MAIN_NAME, step_id]
        process = translator.gen_process_from_cmdtool(tool, sources, scope)
        actual_outputs = {out.get_string() for out in process.outputs}
        expected_outputs = {'path params.stp4_output_filename, emit: out'}
        self.assertEqual(actual_outputs, expected_outputs)
        
    def test_input_selector_array(self) -> None:
        wf = InputSelectorTestWF()
        refresh_workflow_inputs(wf)
        tool = wf.step_nodes["stp3"].tool
        sources = wf.step_nodes["stp3"].sources
        process = translator.gen_process_from_cmdtool(tool, sources, scope=['stp3'])
        actual_outputs = {out.get_string() for out in process.outputs}
        expected_outputs = {'path inp, emit: out'}
        self.assertEqual(actual_outputs, expected_outputs)
    
    def test_input_selector_filename_reference(self) -> None:
        step_id = "stp3"
        tool = self.wf.step_nodes[step_id].tool
        sources = self.wf.step_nodes[step_id].sources
        scope = [settings.NF_MAIN_NAME, step_id]
        process = translator.gen_process_from_cmdtool(tool, sources, scope)
        actual_outputs = {out.get_string() for out in process.outputs}
        expected_outputs = {'path "${inp.simpleName}.recalibrated.bam", emit: out'}
        self.assertEqual(actual_outputs, expected_outputs)
    
    def test_input_selector_filename_generated(self) -> None:
        step_id = "stp2"
        tool = self.wf.step_nodes[step_id].tool
        sources = self.wf.step_nodes[step_id].sources
        scope = [settings.NF_MAIN_NAME, step_id]
        process = translator.gen_process_from_cmdtool(tool, sources, scope)
        actual_outputs = {out.get_string() for out in process.outputs}
        expected_outputs = {'path "generated.recalibrated.bam", emit: out'}
        self.assertEqual(actual_outputs, expected_outputs)
    
    def test_file_pair(self) -> None:
        # eg read1.fastq, read2.fastq
        # collection method is list, len(list) == 2.
        step_id = "stp6"
        tool = self.wf.step_nodes[step_id].tool
        sources = self.wf.step_nodes[step_id].sources
        scope = [settings.NF_MAIN_NAME, step_id]
        process = translator.gen_process_from_cmdtool(tool, sources, scope)
        actual_outputs = {out.get_string() for out in process.outputs}
        expected_outputs = {'path "[${inp.simpleName + \'-R1.fastq\'}, ${inp.simpleName + \'-R2.fastq\'}]", emit: out'}
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
        
    def test_complex_expression(self) -> None:
        # two_value operator etc. uses ${} syntax around whole phrase.
        # strings inside are quoted. 
        step_id = "stp5"
        tool = self.wf.step_nodes[step_id].tool
        sources = self.wf.step_nodes[step_id].sources
        scope = [settings.NF_MAIN_NAME, step_id]
        process = translator.gen_process_from_cmdtool(tool, sources, scope)
        actual_outputs = {out.get_string() for out in process.outputs}
        expected_outputs = {'path "${inp.simpleName + \'.gz\'}", emit: out'}
        self.assertEqual(actual_outputs, expected_outputs)
    
    def test_pythontool(self) -> None:
        wf = OutputsPythonToolTestWF()
        refresh_workflow_inputs(wf)
        
        # file output
        tool = wf.step_nodes["stp0"].tool
        sources = wf.step_nodes["stp0"].sources
        process = translator.gen_process_from_codetool(tool, sources, scope=['stp0'])
        actual_outputs = {out.get_string() for out in process.outputs}
        expected_outputs = {'path "out_out", emit: out'}
        self.assertEqual(actual_outputs, expected_outputs)
        
        # String output
        tool = wf.step_nodes["stp1"].tool
        sources = wf.step_nodes["stp1"].sources
        process = translator.gen_process_from_codetool(tool, sources, scope=['stp1'])
        actual_outputs = {out.get_string() for out in process.outputs}
        expected_outputs = {'val "${file("${task.workDir}/out_out").text}", emit: out'}
        self.assertEqual(actual_outputs, expected_outputs)
        
        # Array(String) output
        tool = wf.step_nodes["stp2"].tool
        sources = wf.step_nodes["stp2"].sources
        process = translator.gen_process_from_codetool(tool, sources, scope=['stp2'])
        actual_outputs = {out.get_string() for out in process.outputs}
        expected_outputs = {'val "${file("${task.workDir}/out_out").text.split(\',\')}", emit: out'}
        self.assertEqual(actual_outputs, expected_outputs)

    @unittest.skip('not implemented')
    def test_translate_commandtool(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_translate_pythontool(self) -> None:
        raise NotImplementedError



class TestProcessScript(unittest.TestCase):
    """
    Tests of the process prescript and script sections.
    """

    def setUp(self) -> None:
        settings.MINIMAL_PROCESS = True
        settings.MODE = 'workflow'

    def test_components_prescript(self) -> None:
        wf = StepInputsTestWF()
        refresh_workflow_inputs(wf)
        tool = wf.step_nodes["stp1"].tool
        sources = wf.step_nodes["stp1"].sources
        process = translator.gen_process_from_cmdtool(tool, sources, scope=['stp1'])
        actual_prescript = process.pre_script
        expected_lines = {
            'def pos_default = params.in_int ? params.in_int : 95',
            'def pos_optional = params.in_str ? params.in_str : ""',
            'def flag_true = params.in_bool == false ? "" : "--flag-true"',
            'def flag_false = params.in_bool ? "--flag-false" : ""',
            'def opt_default = params.in_int ? params.in_int : 5',
            'def opt_optional = params.in_str ? "--opt-optional ${params.in_str}" : ""',
        }
        for ln in expected_lines:
            self.assertIn(ln, actual_prescript)
    
    def test_components_script(self) -> None:
        wf = StepInputsTestWF()
        refresh_workflow_inputs(wf)
        tool = wf.step_nodes["stp1"].tool
        sources = wf.step_nodes["stp1"].sources
        process = translator.gen_process_from_cmdtool(tool, sources, scope=['stp1'])
        actual_script = process.script
        expected_lines = [
            'echo',
            '${pos_basic}',
            '${pos_default}',
            '${pos_optional}',
            '${flag_true}',
            '${flag_false}',
            '--opt-basic ${params.in_str}',
            '--opt-default ${opt_default}',
            '${opt_optional}',
        ]
        for ln in expected_lines:
            self.assertIn(ln, actual_script)

    def test_components_array_prescript(self) -> None:
        wf = ArrayStepInputsTestWF()
        refresh_workflow_inputs(wf)
        tool = wf.step_nodes["stp1"].tool
        sources = wf.step_nodes["stp1"].sources
        process = translator.gen_process_from_cmdtool(tool, sources, scope=['stp1'])
        actual_prescript = process.pre_script
        expected_lines = [
            'def pos_basic = pos_basic.join(\' \')',
            'def pos_basic2 = pos_basic2 ? pos_basic2.join(\' \') : ""',
            'def pos_default = params.in_int_array ? params.in_int_array.join(\' \') : "1 2 3"',
            'def pos_optional = params.in_str_array ? params.in_str_array.join(\' \') : ""',
            'def opt_basic = params.in_str_array.join(\' \')',
            'def opt_default = params.in_int_array ? params.in_int_array.collect{ "--opt-default " + it }.join(\' \') : "--opt-default 1 --opt-default 2 --opt-default 3"',
            'def opt_optional = params.in_str_array ? "--opt-optional," + params.in_str_array.join(\',\') : ""',
            
        ]
        print(actual_prescript)
        for ln in expected_lines:
            self.assertIn(ln, actual_prescript)
    
    def test_components_array_script(self) -> None:
        wf = ArrayStepInputsTestWF()
        refresh_workflow_inputs(wf)
        tool = wf.step_nodes["stp1"].tool
        sources = wf.step_nodes["stp1"].sources
        process = translator.gen_process_from_cmdtool(tool, sources, scope=['stp1'])
        actual_script = process.script
        expected_lines = [
            'echo',
            '${pos_basic}',
            '${pos_basic2}',
            '${pos_default}',
            '${pos_optional}',
            '--opt-basic=${opt_basic}',
            '--opt-default ${opt_default}',
            '${opt_optional}',
        ]
        print(actual_script)
        for ln in expected_lines:
            self.assertIn(ln, actual_script)
    
    def test_secondaries(self) -> None:
        # name accession should be different?
        wf = SecondariesIOTestWF()
        refresh_workflow_inputs(wf)
        tool = wf.step_nodes["stp1"].tool
        sources = wf.step_nodes["stp1"].sources
        process = translator.gen_process_from_cmdtool(tool, sources, scope=['stp1'])
        
        # pre-script
        actual_pre_script = process.pre_script
        expected_pre_script = {}
        for ln in expected_pre_script:
            self.assertIn(ln, actual_pre_script)
        
        # script
        actual_script = process.script
        expected_script = {
            'echo',
            '${bam}',
        }
        for ln in expected_script:
            self.assertIn(ln, actual_script)
    
    @unittest.skip('not implemented')
    def test_secondaries_array(self) -> None:
        # unnecessary as format follows non-array secondaries (test above)
        raise NotImplementedError

    def test_filename_generated_tool(self):
        wf = FilenameGeneratedTestWF()
        refresh_workflow_inputs(wf)
        step_id = 'stp1'
        tool = wf.step_nodes[step_id].tool
        sources = wf.step_nodes[step_id].sources
        scope: list[str] = [settings.NF_MAIN_NAME, step_id]
        process = translator.gen_process_from_cmdtool(tool, sources, scope)
        print(process.get_string())
        expected = """\
process STP1 {
    publishDir "${params.outdir}/stp1"
    debug true

    input:
    path fileInp
    path fileInpOptional

    output:
    val "*", emit: out

    script:
    \"\"\"
    echo \\
    ${params.in_str} \\
    ${params.in_str_opt} \\
    ${fileInp.simpleName}.transformed.fnp \\
    ${fileInpOptional.simpleName}.optional.txt \\
    \"\"\"

}
"""
        self.assertEqual(expected, process.get_string())

    def test_pythontool(self) -> None:
        wf = InputsPythonToolTestWF()
        refresh_workflow_inputs(wf)
        
        # File, String, Int input types
        tool = wf.step_nodes["stp0"].tool
        sources = wf.step_nodes["stp0"].sources
        process = translator.gen_process_from_codetool(tool, sources, scope=['stp0'])
        actual_script = process.script
        expected_lines = {
            'result = code_block(inp1="${inp1}", inp2="${params.in_str}", inp3=${params.in_int})',
        }
        print(process.get_string())
        for ln in expected_lines:
            self.assertIn(ln, actual_script)
        
        # Array(String) input type
        tool = wf.step_nodes["stp1"].tool
        sources = wf.step_nodes["stp1"].sources
        process = translator.gen_process_from_codetool(tool, sources, scope=['stp1'])
        actual_script = process.script
        expected_lines = {
            'result = code_block(inp="${params.in_str_arr}".split(" "))',
        }
        print(process.get_string())
        for ln in expected_lines:
            self.assertIn(ln, actual_script)
        
        # File (secondaries) input type
        tool = wf.step_nodes["stp2"].tool
        sources = wf.step_nodes["stp2"].sources
        process = translator.gen_process_from_codetool(tool, sources, scope=['stp2'])
        actual_script = process.script
        expected_lines = {
            'result = code_block(inp="${bam}")',
        }
        print(process.get_string())
        for ln in expected_lines:
            self.assertIn(ln, actual_script)

    def test_filename_types(self) -> None:
        wf = FilenameTestWF()
        refresh_workflow_inputs(wf)

        step_id = 'stp1'
        tool = wf.step_nodes[step_id].tool
        sources = wf.step_nodes[step_id].sources
        process = translator.gen_process_from_cmdtool(tool, sources, scope=[settings.NF_MAIN_NAME, step_id])
        actual_script = process.script
        expected_lines = {
            'echo',
            '${inp1}',
            '${inp2}',
        }
        print(process.get_string())
        for ln in expected_lines:
            self.assertIn(ln, actual_script)

        step_id = 'stp2'
        tool = wf.step_nodes[step_id].tool
        sources = wf.step_nodes[step_id].sources
        process = translator.gen_process_from_cmdtool(tool, sources, scope=[settings.NF_MAIN_NAME, step_id])
        actual_script = process.script
        expected_lines = {
            'echo',
            '${inp1}',
            '${inp1.simpleName}.processed.txt',
        }
        print(process.get_string())
        for ln in expected_lines:
            self.assertIn(ln, actual_script)

    @unittest.skip('not implemented')
    def test_translate_commandtool(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_translate_pythontool(self) -> None:
        raise NotImplementedError




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
        step_id = 'stp1'
        step = wf.step_nodes[step_id]
        scope = [settings.NF_MAIN_NAME, step_id]
        actual = nfgen.get_args(step, scope)
        expected = [
            "ch_in_file",
            "ch_in_file_opt",
        ]
        self.assertEqual(expected, actual)
        
    # static step inputs
    def test_static_inputs(self):
        wf = StepInputsTestWF()
        refresh_workflow_inputs(wf)
        step_id = 'stp1'
        step = wf.step_nodes[step_id]
        scope = [settings.NF_MAIN_NAME, step_id]
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
        actual = nfgen.get_args(step, scope)
        for tinput_name in not_expected:
            self.assertNotIn(tinput_name, actual)

    # connections
    def test_connections_files(self) -> None:
        wf = StepConnectionsTestWF()
        refresh_workflow_inputs(wf)
        step_id = 'stp2'
        step = wf.step_nodes[step_id]
        scope = [settings.NF_MAIN_NAME, step_id]
        expected = ["STP1.out.out"]
        actual = nfgen.get_args(step, scope)
        self.assertEqual(expected, actual)

    def test_filename_types(self) -> None:
        wf = FilenameTestWF()
        refresh_workflow_inputs(wf)
        step_id = 'stp1'
        step = wf.step_nodes[step_id]
        scope = [settings.NF_MAIN_NAME, step_id]
        actual = nfgen.get_args(step, scope)
        expected = [
            "ch_in_file",
            "ch_in_str",
        ]
        self.assertEqual(expected, actual)
        
        step_id = 'stp2'
        step = wf.step_nodes[step_id]
        scope = [settings.NF_MAIN_NAME, step_id]
        actual = nfgen.get_args(step, scope)
        expected = [
            "ch_in_file",
        ]
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
        step_id = 'stp2'
        step = wf.step_nodes[step_id]
        scope = [settings.NF_MAIN_NAME, step_id]
        actual = nfgen.get_args(step, scope)
        expected = [
            "STP1.out.out"
        ]
        self.assertEqual(expected, actual)
    
    def test_workflow_inputs_array(self) -> None:
        wf = ArrayStepInputsTestWF()
        refresh_workflow_inputs(wf)
        step_id = 'stp1'
        step = wf.step_nodes[step_id]
        scope = [settings.NF_MAIN_NAME, step_id]
        actual = nfgen.get_args(step, scope)
        expected = [
            "ch_in_file_array",
            "ch_in_file_array_opt",
        ]
        self.assertEqual(expected, actual)

    def test_static_step_inputs_array(self):
        wf = ArrayStepInputsTestWF()
        refresh_workflow_inputs(wf)
        step_id = 'stp2'
        step = wf.step_nodes[step_id]
        scope = [settings.NF_MAIN_NAME, step_id]
        actual = nfgen.get_args(step, scope)
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
        step_id = 'stp1'
        step = wf.step_nodes[step_id]
        scope = [settings.NF_MAIN_NAME, step_id]
        actual = nfgen.get_args(step, scope)
        expected = [
            "ch_in_file_array.flatten()"
        ]
        self.assertEqual(expected, actual)
    
    def test_scatter_dot(self) -> None:
        wf = ScatterDotTestWF()
        refresh_workflow_inputs(wf)
        step_id = 'stp1'
        step = wf.step_nodes[step_id]
        scope = [settings.NF_MAIN_NAME, step_id]
        actual = nfgen.get_args(step, scope)
        expected = [
            "ch_in_file_array.flatten()",
            "ch_in_str_array.flatten()",
        ]
        self.assertEqual(expected, actual)
    
    def test_scatter_cross(self) -> None:
        wf = ScatterCrossTestWF()
        refresh_workflow_inputs(wf)
        step_id = 'stp1'
        step = wf.step_nodes[step_id]
        scope = [settings.NF_MAIN_NAME, step_id]
        actual = nfgen.get_args(step, scope)

        # cartesian cross channel manipulation in workflow
        operation = nfgen.channels.gen_scatter_cross_operation(step.sources, step.scatter)
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
        expected = [
            "ch_cartesian_cross.in_file_array",
            "ch_cartesian_cross.in_str_array",
        ]
        self.assertEqual(expected, actual)
    
    def test_scatter_connection(self) -> None:
        wf = ChainedScatterTestWF()
        refresh_workflow_inputs(wf)
        # single -> single
        step_id = 'stp2'
        step = wf.step_nodes[step_id]
        scope = [settings.NF_MAIN_NAME, step_id]
        actual = nfgen.get_args(step, scope)
        expected = [
            "STP1.out.out"
        ]
        self.assertEqual(expected, actual)
        # array -> single
        step_id = 'stp4'
        step = wf.step_nodes[step_id]
        scope = [settings.NF_MAIN_NAME, step_id]
        actual = nfgen.get_args(step, scope)
        expected = [
            "STP3.out.out.flatten()"
        ]
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
        inp = wf.input_nodes['inAlignments']
        expected_channel = 'ch_in_alignments = Channel.fromPath( params.in_alignments_bam, params.in_alignments_bai ).collect()'
        channel = nfgen.channels.getall(inp.uuid)[0]
        actual_channel = channel.declaration
        self.assertEquals(actual_channel, expected_channel)
        
        # step inputs created correctly?
        step_id = 'stp1'
        step = wf.step_nodes[step_id]
        scope = [settings.NF_MAIN_NAME, step_id]
        actual_inputs = nfgen.get_args(step, scope)
        expected_inputs = [
            "ch_in_alignments"
        ]
        self.assertEquals(actual_inputs, expected_inputs)
    
    def test_secondaries_connections(self) -> None:
        wf = SecondariesConnectionsTestWF()
        refresh_workflow_inputs(wf)
        # step inputs created correctly?
        step_id = 'stp2'
        step = wf.step_nodes[step_id]
        scope = [settings.NF_MAIN_NAME, step_id]
        actual_inputs = nfgen.get_args(step, scope)
        expected_inputs = [
            "STP1.out.out"
        ]
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
    # 'in_bool_array': fr'\[\]',
    # 'stp2_flag_true': fr'\[true\]',
    # 'stp2_flag_false': fr'\[true\]',
    'in_file_array': fr"""\[
        // list files here
    \]""",
    'in_file_array_opt': fr"""\[
        // list files here
    \]""",
    
    'in_str_array': fr"""\[
        // list strings here
    \]""",
    'in_int_array': fr'\[\]',

    'stp2_pos_default': fr'\[4, 5, 6\]',
    'stp2_pos_optional': fr"""\[
        'hi',
        'there',
        'friend',
    \]""",
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
        wf = ConditionStepTestWF()
        refresh_workflow_inputs(wf)
        step_id = 'print'
        step = wf.step_nodes[step_id]
        scope = [settings.NF_MAIN_NAME, step_id]
        actual = nfgen.get_args(step, scope)
        expected = [
            "[params.mystring, GET_STRING.out.out].first()"
        ]
        self.assertEqual(actual, expected)

    def test_with_expression(self):
        wf = StepInputExpressionTestWF()
        refresh_workflow_inputs(wf)
        step_id = 'print'
        step = wf.step_nodes[step_id]
        scope = [settings.NF_MAIN_NAME, step_id]
        actual = nfgen.get_args(step, scope)
        expected = [
            "binding.hasVariable(params.mystring) ? params.mystring : params.mystring_backup"
        ]
        self.assertEqual(actual, expected)



class TestUnwrap(unittest.TestCase):

    def setUp(self) -> None:
        settings.MINIMAL_PROCESS = True
        settings.MODE = 'workflow'
        wf = UnwrapTestWF()
        refresh_workflow_inputs(wf)
        step_id = 'stp1'
        tool = wf.step_nodes[step_id].tool
        sources = wf.step_nodes[step_id].sources
        self.prescript, self.script = nfgen.process.gen_script_for_cmdtool(
            tool=tool,
            scope=[settings.NF_MAIN_NAME, step_id],
            sources=sources,
            stdout_filename='out'
        )
        self.wf = StepConnectionsTestWF()
        refresh_workflow_inputs(self.wf)

    def test_filename_generated(self) -> None:
        self.assertIn("--filenameGen myfile.gz", self.script)
    
    @unittest.skip('not implemented')
    def test_filename_input_selector(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_two_value_operator(self) -> None:
        raise NotImplementedError

    def test_input_selector_process_input(self) -> None:
        wf = InputSelectorTestWF()
        refresh_workflow_inputs(wf)
        step_id = 'stp1'
        output_id = 'out'
        process_inputs = {'inp'}
        param_inputs = {}
        internal_inputs = {}
        tool = wf.step_nodes[step_id].tool
        sel = [x.selector for x in tool.outputs() if x.id() == output_id][0]
        actual = nfgen.unwrap_expression(
            val=sel, 
            tool=tool, 
            process_inputs=process_inputs,
            param_inputs=param_inputs,
            internal_inputs=internal_inputs,
            add_curly_braces=True
        )
        expected = '${inp}'
        self.assertEqual(actual, expected)

    def test_input_selector_param_input(self) -> None:
        wf = InputSelectorTestWF()
        refresh_workflow_inputs(wf)
        step_id = 'stp2'
        output_id = 'out'
        process_inputs = {}
        param_inputs = {'inp'}
        internal_inputs = {}
        tool = wf.step_nodes[step_id].tool
        sources = wf.step_nodes[step_id].sources
        sel = [x.selector for x in tool.outputs() if x.id() == output_id][0]
        actual = nfgen.unwrap_expression(
            val=sel, 
            tool=tool, 
            sources=sources,
            process_inputs=process_inputs,
            param_inputs=param_inputs,
            internal_inputs=internal_inputs,
            add_curly_braces=True
        )
        expected = "${params.in_str}"
        self.assertEqual(actual, expected)
    
    def test_input_node_channel(self) -> None:
        # file input (channel)
        node = self.wf.input_nodes['inFile']
        actual = nfgen.unwrap_expression(node)
        expected = 'ch_in_file'
        self.assertEqual(actual, expected)

    def test_input_node_param(self) -> None:
        # nonfile input (param)
        node = self.wf.input_nodes['inStr']
        actual = nfgen.unwrap_expression(node)
        expected = 'params.in_str'
        self.assertEqual(actual, expected)

    def test_step_connection(self) -> None:
        wf = StepConnectionsTestWF()
        refresh_workflow_inputs(wf)
        step_id = "stp2"
        inp_id = 'inp'
        sources = self.wf.step_nodes[step_id].sources
        src = sources[inp_id]
        actual = nfgen.unwrap_expression(src)
        expected = 'STP1.out.out'
        self.assertEqual(actual, expected)
    
    def test_alias_selector(self) -> None:
        wf = AliasSelectorTestWF()
        refresh_workflow_inputs(wf)
        step_id = "stp2"
        inp_id = 'inp'
        sources = wf.step_nodes[step_id].sources
        src = sources[inp_id]
        actual = nfgen.unwrap_expression(src)
        expected = 'STP1.out.out'
        self.assertEqual(actual, expected)
    
    def test_first_operator(self) -> None:
        wf = ConditionStepTestWF()
        refresh_workflow_inputs(wf)
        step_id = "print"
        inp_id = 'inp'
        sources = wf.step_nodes[step_id].sources
        src = sources[inp_id]
        actual = nfgen.unwrap_expression(src)
        expected = '[params.mystring, GET_STRING.out.out].first()'
        self.assertEqual(actual, expected)
    
    def test_index_operator(self) -> None:
        wf = IndexOperatorTestWF()
        refresh_workflow_inputs(wf)
        step_id = "stp1"
        inp_id = 'inp'
        sources = wf.step_nodes[step_id].sources
        src = sources[inp_id]
        actual = nfgen.unwrap_expression(src)
        expected = 'ch_in_file_arr[0]'
        self.assertEqual(actual, expected)




class TestStringFormatter(unittest.TestCase):
    
    def setUp(self) -> None:
        settings.MINIMAL_PROCESS = True
        settings.MODE = 'workflow'
    
    def test_string_formatter(self):
        tool = BasicTestTool()
        sf = StringFormatter("no format")
        res = nfgen.unwrap_expression(sf, tool)
        self.assertEqual("no format", res)

    def test_string_formatter_string(self):
        tool = BasicTestTool()
        sf = StringFormatter("there's a {str_arg} arg", str_arg="string")
        res = nfgen.unwrap_expression(sf, tool)
        self.assertEqual("there's a ${\"string\"} arg", res)
    
    def test_string_formatter_inputselector_process_input(self):
        tool = BasicTestTool()
        sf = StringFormatter("an input {arg}", arg=InputSelector("testtool"))
        process_inputs = {'testtool'}
        param_inputs = {}
        internal_inputs = {}
        actual = nfgen.unwrap_expression(
            val=sf, 
            tool=tool,
            process_inputs=process_inputs,
            param_inputs=param_inputs,
            internal_inputs=internal_inputs,
        )
        expected = 'an input ${testtool}'
        self.assertEqual(actual, expected)
    
    def test_string_formatter_inputselector_param_input(self):
        wf = StringFormatterTestWF()
        refresh_workflow_inputs(wf)
        step_id = 'stp1'
        tool = wf.step_nodes[step_id].tool
        sources = wf.step_nodes[step_id].sources
        process_inputs = {}
        param_inputs = {'testtool'}
        internal_inputs = {}
        sf = StringFormatter("an input {arg}", arg=InputSelector("testtool"))
        actual = nfgen.unwrap_expression(
            val=sf, 
            tool=tool,
            sources=sources,
            process_inputs=process_inputs,
            param_inputs=param_inputs,
            internal_inputs=internal_inputs,
        )
        expected = 'an input ${params.in_str}'
        self.assertEqual(actual, expected)

    def test_string_formatter_two_param(self):
        tool = InputQualityTestTool()
        sf = StringFormatter(
            "{username}:{password}",
            username=InputSelector("user"),
            password=InputSelector("static"),
        )
        sources = {}
        process_inputs = {'user', 'static'}
        param_inputs = {}
        internal_inputs = {}
        actual = nfgen.unwrap_expression(
            val=sf,
            tool=tool,
            sources=sources,
            process_inputs=process_inputs,
            param_inputs=param_inputs,
            internal_inputs=internal_inputs,
        )
        expected = '${user}:${static}'
        self.assertEqual(actual, expected)

    def test_escaped_characters(self):
        tool = InputQualityTestTool()
        sf = StringFormatter(
            "{username}\\t{password}",
            username=InputSelector("user"),
            password=InputSelector("static"),
        )
        sources = {}
        process_inputs = {'user', 'static'}
        param_inputs = {}
        internal_inputs = {}
        actual = nfgen.unwrap_expression(
            val=sf,
            tool=tool,
            sources=sources,
            process_inputs=process_inputs,
            param_inputs=param_inputs,
            internal_inputs=internal_inputs,
        )
        expected = '${user}\\t${static}'
        self.assertEqual(actual, expected)

    def test_expression_arg(self):
        tool = BasicTestTool()
        sf = StringFormatter(
            "{name}:{items}",
            name=InputSelector("testtool"),
            items=JoinOperator(InputSelector("arrayInp"), separator=";"),
        )
        sources = {}
        process_inputs = {'testtool', 'arrayInp'}
        param_inputs = {}
        internal_inputs = {}
        actual = nfgen.unwrap_expression(
            val=sf,
            tool=tool,
            sources=sources,
            process_inputs=process_inputs,
            param_inputs=param_inputs,
            internal_inputs=internal_inputs,
        )
        expected = "${testtool}:${arrayInp.join(\";\")}"
        self.assertEqual(actual, expected)



class TestSubWorkflows(unittest.TestCase):
    # sometimes logic here a bit more complex or weird 
    # due to subworkflows. apologies. 
    
    def setUp(self) -> None:
        settings.MINIMAL_PROCESS = True
        settings.MODE = 'workflow'
        self.wf = SubworkflowTestWF()
        refresh_workflow_inputs(self.wf)

    @unittest.skip('not implemented')
    def test_param_system(self) -> None:
        # currently params system doesnt reach to subworkflows. 
        # will implement in future (time permitting)
        raise NotImplementedError
    
    def test_channel_system(self) -> None:
        step_id = 'apples_subworkflow'
        subwf = self.wf.step_nodes[step_id].tool
        for inp in subwf.input_nodes.values():
            relevant_channel = nfgen.channels.get(inp.uuid)
            assert(relevant_channel)   # 1 channel per each subworkflow input
    
    def test_files_created(self) -> None:
        mainstr, substr_dict = translator.translate_workflow(self.wf)
        expected_filepaths = set([
            'modules/file_tool',
            'modules/string_tool',
            'modules/int_tool',
            'subworkflows/oranges_subworkflow',
            'subworkflows/apples_subworkflow',
        ])
        actual_filepaths = set(substr_dict.keys())
        self.assertEqual(actual_filepaths, expected_filepaths)

    def test_structure(self) -> None:
        # take main emit
        mainstr, substr_dict = translator.translate_workflow(self.wf)
        self.assertNotIn('take:', mainstr)
        self.assertNotIn('main:', mainstr)
        self.assertNotIn('emit:', mainstr)
        subwfstr = substr_dict['subworkflows/oranges_subworkflow']
        self.assertIn('take:', subwfstr)
        self.assertIn('main:', subwfstr)
        self.assertIn('emit:', subwfstr)

    def test_call(self) -> None:
        # translate workflow, building all nf items and files
        translator.translate_workflow(self.wf)

        # focusing in on specific subworkflow
        step_id = 'apples_subworkflow'
        step = self.wf.step_nodes[step_id]
        scope = [settings.NF_MAIN_NAME, step_id]
        args = nfgen.get_args(step, scope)

        # generate nf subworkflow object
        nf_workflow = translator.gen_workflow(
            name=step_id, 
            scope=scope,
            sources=step.sources,
            wf=step.tool
        )
        
        # check the arg order matches the subworkflow input channel order
        expected_arg_order = ['params.in_int', 'params.in_str']
        expected_channel_order = ['ch_in_int', 'ch_in_str']
        actual_arg_order = args
        actual_channel_order = [t.get_string() for t in nf_workflow.take]

        self.assertEquals(actual_arg_order, expected_arg_order)
        self.assertEquals(actual_channel_order, expected_channel_order)

    def test_imports(self) -> None:
        # translate workflow, building all nf items and files
        wf = SubworkflowTestWF()
        translator.translate_workflow(wf)

        # focusing in on specific subworkflow
        step_id = 'apples_subworkflow'
        scope = [settings.NF_MAIN_NAME, step_id]
        subwf_file = translator.file_register.get(scope)

        self.assertEquals(len(subwf_file.imports), 2)

        expected_imports = [
            {
                'name': 'string_tool', 
                'source': 'modules/string_tool'
            },
            {
                'name': 'oranges_subworkflow', 
                'source': 'subworkflows/oranges_subworkflow'
            },
        ]
        actual_imports = []
        for nf_import in subwf_file.imports:
            import_item = {
                'name': nf_import.items[0].name, 
                'source': nf_import.source
            },
            actual_imports += import_item
        
        self.assertEquals(actual_imports, expected_imports)
    
    @unittest.skip('not implemented')
    def test_nested(self) -> None:
        raise NotImplementedError
    
