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
    SecondaryInputTestTool,
    BasicTestTool,
    AppendedSecondaryOutputTestTool,
    ReplacedSecondaryOutputTestTool,
)

from janis_core.tests.testworkflows import (

    # basics
    BasicIOTestWF,
    StepInputsTestWF,
    StepInputsWFInputTestWF,
    StepInputsStaticTestWF,
    StepInputsPartialStaticTestWF,
    StepInputsMinimalTestWF,
    StepConnectionsTestWF,

    # arrays
    ArrayIOTestWF,
    ArrayStepInputsTestWF,
    ArrayStepConnectionsTestWF,

    # scatter
    BasicScatterTestWF,
    ChainedScatterTestWF,
    MultiFieldScatterTestWF,

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
from janis_core import Array, String, File, Boolean, Filename



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

def register_params(wf: Workflow) -> dict[str, str]:
    scope: list[str] = []
    nfgen.params.clear()
    nfgen.params.register(the_entity=wf, scope=scope)
    for step in wf.step_nodes.values():
        current_scope = deepcopy(scope)
        current_scope.append(step.id())
        nfgen.params.register(the_entity=step.tool, sources=step.sources, scope=current_scope)
    params = nfgen.params.getall()
    return {p.name: p.value for p in params} 

def register_params_tool(tool: CommandTool) -> dict[str, str]:
    scope: list[str] = []
    nfgen.params.clear()
    nfgen.params.register(the_entity=tool, scope=scope)
    params = nfgen.params.getall()
    return {p.name: p.value for p in params} 

def register_channels(wf: Workflow) -> list[Channel]:
    nfgen.channels.clear()
    nfgen.channels.register(wf)
    return nfgen.channels.getall()



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
        register_params(wf)
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
        settings.MINIMAL_PROCESS = False
        actual_params = register_params(wf)
        expected_params = {
            'in_file': 'null',
            'stp1_pos_default': 'null',
            'stp1_pos_optional': 'null',
            'stp1_flag_true': 'null',
            'stp1_flag_false': 'null',
            'stp1_opt_basic': 'null',
            'stp1_opt_default': 'null',
            'stp1_opt_optional': 'null'
        }
        self.assertEquals(actual_params, expected_params)
        
        # minimal
        settings.MINIMAL_PROCESS = True
        actual_params = register_params(wf)
        expected_params = {
            'in_file': 'null',
            'stp1_pos_default': 'null',
            'stp1_pos_optional': 'null',
            'stp1_flag_true': 'null',
            'stp1_flag_false': 'null',
            'stp1_opt_basic': 'null',
            'stp1_opt_default': 'null',
            'stp1_opt_optional': 'null'
        }
        self.assertEquals(actual_params, expected_params)

    def test_static_step_inputs(self) -> None:
        wf = StepInputsStaticTestWF()
        
        # full
        settings.MINIMAL_PROCESS = False
        actual_params = register_params(wf)
        expected_params = {
            'in_file': 'null',
            'stp2_pos_default': '100',
            'stp2_pos_optional': "'static'",
            'stp2_flag_true': 'false',
            'stp2_flag_false': 'true',
            'stp2_opt_basic': "'static'",
            'stp2_opt_default': '100',
            'stp2_opt_optional': "''"
        }
        self.assertEquals(actual_params, expected_params)
        
        # minimal
        settings.MINIMAL_PROCESS = True
        actual_params = register_params(wf)
        expected_params = {
            'in_file': 'null',
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
        
        # full
        settings.MINIMAL_PROCESS = False
        actual_params = register_params(wf)
        expected_params = {
            'in_file': 'null',
            'stp3_pos_default': 'null',
            'stp3_pos_optional': 'null',
            'stp3_flag_true': 'null',
            'stp3_flag_false': 'null',
            'stp3_opt_basic': "'static'",
            'stp3_opt_default': '100',
            'stp3_opt_optional': "''"
        }
        self.assertEquals(actual_params, expected_params)
        
        # minimal
        settings.MINIMAL_PROCESS = True
        actual_params = register_params(wf)
        expected_params = {
            'in_file': 'null',
            'stp3_opt_basic': "'static'",
            'stp3_opt_default': '100',
            'stp3_opt_optional': "''"
        }
        self.assertEquals(actual_params, expected_params)
    
    def test_omitted_step_inputs(self) -> None:
        wf = StepInputsMinimalTestWF()
        
        # full
        settings.MINIMAL_PROCESS = False
        actual_params = register_params(wf)
        expected_params = {
            'in_file': 'null',
            'stp4_pos_default': 'null',
            'stp4_pos_optional': 'null',
            'stp4_flag_true': 'null',
            'stp4_flag_false': 'null',
            'stp4_opt_basic': 'null',
            'stp4_opt_default': 'null',
            'stp4_opt_optional': 'null'
        }
        self.assertEquals(actual_params, expected_params)
        
        # minimal
        settings.MINIMAL_PROCESS = True
        actual_params = register_params(wf)
        expected_params = {
            'in_file': 'null',
            'stp4_opt_basic': 'null',
        }
        self.assertEquals(actual_params, expected_params)

    def test_array_inputs(self) -> None:
        wf = ArrayIOTestWF()
        actual_params = register_params(wf)
        expected_params = {
            'stp1_ins': '[]',
            'in_file_array': '[]',
            'stp3_ins': '[]',
        }
        self.assertEquals(actual_params, expected_params)

    def test_secondaries_inputs(self) -> None:
        wf = SecondariesIOTestWF()
        actual_params = register_params(wf)
        expected_params = {
            'in_alignments_bam': 'null',
            'in_alignments_bai': 'null',
        }
        self.assertEquals(actual_params, expected_params)
    
    def test_secondaries_array_inputs(self) -> None:
        wf = ArraySecondariesTestWF()
        actual_params = register_params(wf)
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
        register_params(wf)
        register_channels(wf)
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
        register_params(wf)
        register_channels(wf)
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
        register_params(wf)
        register_channels(wf)
        nonfile_wf_input_ids = {
            'unicycler_kmers',
            'unicycler_scores',
            'unicycler_startGeneCov',
            'unicycler_startGeneId',
        }
        channel_janis_references = {c.reference for c in nfgen.channels.getall()}
        for winp_id in nonfile_wf_input_ids:
            self.assertNotIn(winp_id, channel_janis_references)

    def test_array_inputs(self) -> None:
        wf = ArrayIOTestWF()
        register_params(wf)
        register_channels(wf)
        channels_ids = {c.name for c in nfgen.channels.getall()}
        expected_ids = {
            'ch_in_file_array',
        }
        self.assertEqual(channels_ids, expected_ids)
        for c in nfgen.channels.getall():
            self.assertTrue(c.collect)
        
    def test_secondaries_inputs(self) -> None:
        wf = SecondariesIOTestWF()
        register_params(wf)
        register_channels(wf)
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
        register_params(wf)
        register_channels(wf)
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

    def test_operator_resource_tool(self):
        # TODO FIX ME
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
    
    @unittest.skip('not implemented')
    def test_metadata(self) -> None:
        # debug true
        # label
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_publishDir(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_container(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_resource_directives(self) -> None:
        # memory, cpus, disk, time
        raise NotImplementedError

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
        from janis_core.tests.testworkflows import AssemblyTestWF
        self.maxDiff = None
        self.wf = AssemblyTestWF()
        register_params(self.wf)

    def test_wf_inputs(self) -> None:
        # need a process input for each File wf input in step sources.
        # non-files are fed data via params.
        step = self.wf.step_nodes['unicycler']
        actual_ids = nfgen.utils.get_process_input_ids(step.tool, step.sources)
        expected_ids = {
            'option1',
            'option2',
            'optionL',
        }
        self.assertEqual(actual_ids, expected_ids)
    
    def test_connections(self) -> None:
        # need a process input for each connection in step sources.
        step = self.wf.step_nodes['CatTestTool']
        actual_ids = nfgen.utils.get_process_input_ids(step.tool, step.sources)
        expected_ids = {
            'inp',
        }
        self.assertEqual(actual_ids, expected_ids)
    
    def test_static_inputs(self) -> None:
        # DO NOT need a process input for each static value in step sources.
        # non-files are fed data via params. 
        step = self.wf.step_nodes['unicycler']
        actual_ids = nfgen.utils.get_process_input_ids(step.tool, step.sources)
        non_expected_ids = {
            'kmers',
            'scores',
            'startGeneCov',
            'startGeneId',
        }
        for input_id in non_expected_ids:
            self.assertNotIn(input_id, actual_ids)
    
    def test_internal_inputs(self) -> None:
        # DO NOT need a process input for each non-exposed tool input. 
        # tool input will be autofilled or ignored in script.
        step = self.wf.step_nodes['unicycler']
        actual_ids = nfgen.utils.get_process_input_ids(step.tool, step.sources)
        expected_ids = {
            'option1',
            'option2',
            'optionL',
        }
        non_expected_ids = {x.id() for x in step.tool.inputs() if x.id() not in expected_ids}
        for input_id in non_expected_ids:
            self.assertNotIn(input_id, actual_ids)

    @unittest.skip('not implemented')
    def test_array_inputs(self) -> None:
        raise NotImplementedError

    @unittest.skip('not implemented')
    def test_secondaries_inputs(self) -> None:
        raise NotImplementedError

    def test_tool_with_secondary_input(self):
        # TODO FIX ME
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
    
    
    @unittest.skip('not implemented')
    def test_secondaries_array_inputs(self) -> None:
        raise NotImplementedError

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
    
    def test_stdout_out_tool(self):
        # TODO FIX ME
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

    def test_tool_with_secondary_output(self):
        # TODO FIX ME
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
        # TODO FIX ME
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
    
    @unittest.skip('not implemented')
    def test_array_outputs(self) -> None:
        raise NotImplementedError

    @unittest.skip('not implemented')
    def test_secondaries_outputs(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_secondaries_array_outputs(self) -> None:
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
    # TODO test more features / edge cases?

    def setUp(self) -> None:
        settings.MINIMAL_PROCESS = True
        settings.MODE = 'workflow'

    def test_basic(self):
        wf = AssemblyTestWF()

        tool = wf.step_nodes["fastqc1"].tool
        sources = wf.step_nodes["fastqc1"].sources
        expected = {
            "inputFile": "ch_in_forward_reads",
            "adapters": "ch_fastqc1_adapters",
            "contaminants": "ch_fastqc1_contaminants",
            "limits": "ch_fastqc1_limits",
        }
        actual = translator.gen_step_inval_dict(tool, sources)
        self.assertEqual(expected, actual)
        
        tool = wf.step_nodes["fastqc3"].tool
        sources = wf.step_nodes["fastqc3"].sources
        expected = {"inputFile": "ch_test_input"}
        actual = translator.gen_step_inval_dict(tool, sources)
        self.assertEqual(expected, actual)
        
        tool = wf.step_nodes["CatTestTool"].tool
        sources = wf.step_nodes["CatTestTool"].sources
        expected = {"inp": "fastqc3.out.outTextFile"}
        actual = translator.gen_step_inval_dict(tool, sources)
        self.assertEqual(expected, actual)
        
        tool = wf.step_nodes["unicycler"].tool
        sources = wf.step_nodes["unicycler"].sources
        expected = {
            "option1": "ch_in_forward_reads",
            "option2": "ch_in_reverse_reads",
            "optionL": "ch_in_long_reads",
        }
        actual = translator.gen_step_inval_dict(tool, sources)
        self.assertEqual(expected, actual)
    
    @unittest.skip('not implemented')
    def test_arrays(self) -> None:
        raise NotImplementedError
    



class TestPlumbingScatter(unittest.TestCase):
    # TODO test more features / edge cases?

    def setUp(self) -> None:
        settings.MINIMAL_PROCESS = True
        settings.MODE = 'workflow'
    
    @unittest.skip('not implemented')
    def test_scatter(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_scatter_dot(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_scatter_cross(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_scatter_array(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_scatter_dot_array(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_scatter_cross_array(self) -> None:
        raise NotImplementedError




class TestPlumbingSecondaries(unittest.TestCase):
    # TODO test more features / edge cases?

    def setUp(self) -> None:
        settings.MINIMAL_PROCESS = True
        settings.MODE = 'workflow'
    
    @unittest.skip('not implemented')
    def test_secondaries_connections(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_secondaries_array_connections(self) -> None:
        raise NotImplementedError
    

    

class TestNextflowConfig(unittest.TestCase):

    def setUp(self) -> None:
        settings.MINIMAL_PROCESS = True
        settings.MODE = 'workflow'
        self.wf = AssemblyTestWF()
        register_params(self.wf)

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
    
    def test_nextflow_config(self) -> None:
        """
        Every wf input should have a param. 
        """
        params = translator.build_inputs_file(self.wf)
        config = translator.stringify_translated_inputs(params)
        expected_params = {
            'in_forward_reads',
            'in_reverse_reads',
            'in_long_reads',
            'test_input',
            'fastqc1_adapters',
            'fastqc1_contaminants',
            'fastqc1_limits',
            'fastqc2_adapters',
            'fastqc2_contaminants',
            'fastqc2_limits',
        }
        for param in expected_params:
            self.assertIn(param, config)

    @unittest.skip('not implemented')
    def test_array_inputs(self) -> None:
        raise NotImplementedError

    @unittest.skip('not implemented')
    def test_secondaries_inputs(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_secondaries_array_inputs(self) -> None:
        raise NotImplementedError

    @unittest.skip('not implemented')
    def test_tool_config(self) -> None:
        raise NotImplementedError

    
    ### OLD

    def test_input_in_input_value_nooptional_nodefault(self):
        wf = WorkflowBuilder(
            "test_nf_input_in_input_value_nooptional_nodefault",
        )
        wf.input("inpId", String())
        ai = {"inpId": "2"}
        self.assertDictEqual(
            {"inpId": "2"}, translator.build_inputs_file(wf, additional_inputs=ai)
        )

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
        )




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
        inputs = translator.gen_step_inval_dict(tool, sources)
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




"""


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