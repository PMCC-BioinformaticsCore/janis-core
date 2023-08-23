
from typing import Any, Optional
import unittest
import os 
import json
import xml.etree.ElementTree as et

from janis_core.ingestion.main import ingest_galaxy

from janis_core.ingestion.galaxy.gxtool.text.simplification.aliases import resolve_aliases
from janis_core.ingestion.galaxy.gxtool.text.simplification.main_statement import mark_main_statement

from janis_core.ingestion.galaxy import runtime
from janis_core.ingestion.galaxy.gxworkflow import load_tool_state
from janis_core.ingestion.galaxy.gxtool.text.cheetah.evaluation import sectional_evaluate
from janis_core.ingestion.galaxy.gxtool.parsing import load_xmltool
from janis_core.ingestion.galaxy.gxtool.text.simplification.simplify import simplify_cmd

from janis_core.ingestion.galaxy.gxworkflow.parsing.tool_step.metadata import parse_step_metadata
from janis_core.ingestion.galaxy.gxtool.model import XMLTool
from janis_core.ingestion.galaxy.gxtool.model import XMLCondaRequirement
from janis_core.ingestion.galaxy.gxtool.command import gen_command

from janis_core.ingestion.galaxy import regex_to_glob
from janis_core.ingestion.galaxy import datatypes
from janis_core.ingestion.galaxy.datatypes.core import file_t, string_t, bool_t

from janis_core import WorkflowBuilder, Workflow
from janis_core import WorkflowMetadata
from janis_core.workflow.workflow import InputNode
from janis_core.workflow.workflow import OutputNode
from janis_core.workflow.workflow import StepNode
from janis_core.workflow.workflow import ScatterMethod

from janis_core import CommandToolBuilder
from janis_core import CommandTool
from janis_core import ToolInput
from janis_core import ToolOutput
from janis_core import ToolMetadata
from janis_core import (
    InputSelector,
    WildcardSelector
)
from janis_core import (
    Stdout,
    Array,
    Float,
    Boolean,
    File,
    String
)

from janis_core import settings
from janis_core.ingestion.galaxy import datatypes
from janis_core.ingestion.galaxy.containers import resolve_dependencies_as_container

from janis_core.ingestion.galaxy.janis_mapping.workflow import to_janis_workflow
from janis_core.ingestion.galaxy.janis_mapping.workflow import to_janis_inputs_dict
from janis_core.ingestion.galaxy.janis_mapping.tool import to_janis_datatype
from janis_core.ingestion.galaxy.janis_mapping.tool import to_janis_selector
from janis_core.ingestion.galaxy.janis_mapping.tool import to_janis_metadata
from janis_core.ingestion.galaxy.janis_mapping.tool import to_janis_tool_input
from janis_core.ingestion.galaxy.janis_mapping.tool import to_janis_tool_output
from janis_core.ingestion.galaxy.janis_mapping.tool import to_janis_tool

# mock objects
from .mock.mock_components import MOCK_POSITIONAL1
from .mock.mock_components import MOCK_FLAG1
from .mock.mock_components import MOCK_OPTION2
from .mock.mock_components import MOCK_REDIRECT_OUTPUT
from .mock.mock_entities import MOCK_WORKFLOW_INPUT1
from .mock.mock_tool import MOCK_TOOL_ABRICATE
from .mock.mock_workflow import MOCK_WORKFLOW
from janis_core.ingestion import ingest
from janis_core.translations import translate

 
QUERY1 = XMLCondaRequirement(_name='abricate', _version='1.0.1')
QUERY1_EXPECTED_RESULT = 'quay.io/biocontainers/abricate:1.0.1--ha8f3691_1'

QUERY2 = XMLCondaRequirement(_name='samtools', _version='1.15')
QUERY2_EXPECTED_RESULT = 'quay.io/biocontainers/samtools:1.15--h1170115_1'

QUERY3 = XMLCondaRequirement(_name='cutadapt', _version='3.5')
QUERY3_EXPECTED_RESULT = 'quay.io/biocontainers/cutadapt:3.5--py36h91eb985_1'

GALAXY_TESTDATA_PATH = os.path.join(os.getcwd(), 'janis_core/tests/data/galaxy')


### helper functions ###

def _reset_global_settings() -> None:
    settings.ingest.galaxy.GEN_IMAGES = False
    settings.ingest.galaxy.DISABLE_CONTAINER_CACHE = False
    settings.ingest.SAFE_MODE = True
    settings.ingest.cwl.INGEST_JAVASCRIPT_EXPRESSIONS = False
    settings.ingest.cwl.REQUIRE_CWL_VERSION = False
    settings.testing.TESTING_USE_DEFAULT_CONTAINER = True
    settings.datatypes.ALLOW_UNPARSEABLE_DATATYPES = True
    settings.graph.ALLOW_UNKNOWN_SOURCE = False
    settings.graph.ALLOW_UNKNOWN_SCATTER_FIELDS = False
    datatypes.populate()

def _run(filepath: str, srcfmt: str, destfmt: str) -> Optional[str]:
    wf = ingest(filepath, srcfmt)
    return translate(wf, destfmt, allow_empty_container=True, export_path='./translated')

def _load_gxworkflow(filepath: str) -> dict[str, Any]:
    with open(filepath, 'r') as fp:
        return json.load(fp)

def _load_tool_state(step: dict[str, Any], additional_filters: list[str]=[]) -> dict[str, Any]:
    metadata = parse_step_metadata(step)
    runtime.tool.update_via_wrapper(metadata.wrapper)
    return load_tool_state(step, additional_filters=additional_filters)

def _configure_tool_settings(step: dict[str, Any]) -> None:
    metadata = parse_step_metadata(step)
    runtime.tool.update_via_wrapper(metadata.wrapper)

def _read_cmd(path: str) -> str:
    tree = et.parse(path)
    root = tree.getroot()
    assert(root.text)
    return root.text

def _load_xmltool_for_step(filepath: str, step: int) -> XMLTool:
    gx_workflow = _load_gxworkflow(filepath)
    gx_step = gx_workflow['steps'][str(step)]
    _configure_tool_settings(gx_step)
    return load_xmltool(runtime.tool.tool_path)

def _docker_not_running() -> bool:
    import subprocess
    try:
        completed_process = subprocess.run(['docker', 'version'], shell=True, capture_output=True)
        if completed_process.returncode == 0:
            return True
        else:
            return False
    except FileNotFoundError:
        return False
    

# def _prepare_tool_state_for_cheetah(path: str, step: int) -> dict[str, Any]:
#     with open(path, 'r') as fp:
#         gxworkflow = json.load(fp)
#     step = gxworkflow['steps'][str(step)]
#     _configure_tool_settings(step)
#     tool_state = load_tool_state(step)
#     return tool_state

# def _load_tool_command_and_state(step: dict[str, Any], flat: bool=False) -> Tuple[str, dict[str, Any]]:
#     xmltool = load_xmltool(runtime.tool.tool_path)
#     tool_state = load_tool_state(step, flat=flat)
#     cmdstr = simplify_xmltool_command(
#         xmltool=xmltool, 
#         inputs_dict=tool_state,
#         # additional_filters=['remove_dataset_attributes', 'remove_dataset_methods']
#     )
#     return cmdstr, tool_state



### test classes ###

class TestRegexToGlob(unittest.TestCase):

    def setUp(self) -> None:
        _reset_global_settings()
    
    def test_convert_wildcards(self) -> None:
        self.assertEqual(regex_to_glob.convert('hello .*?'), 'hello *')
        self.assertEqual(regex_to_glob.convert('my [^_]+ string'), 'my * string')
        self.assertEqual(regex_to_glob.convert('(0|1)+asd'), '*asd')
        self.assertEqual(regex_to_glob.convert(' \S+ \w+? single'), ' * * single')
    
    def test_convert_logical_or(self) -> None:
        self.assertEqual(regex_to_glob.convert('(.*\.tar|.*\.gz)'), '{*.tar,*.gz}')
        self.assertEqual(regex_to_glob.convert(' (cat|dog|bat) '), ' {cat,dog,bat} ')
    
    def test_convert_char_set(self) -> None:
        self.assertEqual(regex_to_glob.convert('[123]'), '[1,2,3]')
        self.assertEqual(regex_to_glob.convert('[a-z123]'), '[a-z,1,2,3]')
        self.assertEqual(regex_to_glob.convert('[^a-m]'), '[!a-m]')
        self.assertEqual(regex_to_glob.convert('[a-z]'), '[a-z]')
    
    def test_convert_special_chars(self) -> None:
        self.assertEqual(regex_to_glob.convert('hello \S there \w '), 'hello ? there ? ')
        self.assertEqual(regex_to_glob.convert('hello . there'), 'hello ? there')

    def test_convert_easy(self) -> None:
        self.assertEqual(regex_to_glob.convert('report_data/multiqc_.+\.txt'), 'report_data/multiqc_*.txt')
        self.assertEqual(regex_to_glob.convert('.+\.png'), '*.png')
        self.assertEqual(regex_to_glob.convert('mqc_.+\.txt'), 'mqc_*.txt')
    
    def test_convert_medium(self) -> None:
        self.assertEqual(regex_to_glob.convert(".+vs_.+)"), "*vs_*)") 
        self.assertEqual(regex_to_glob.convert(".+_ss\.ps"), "*_ss.ps") 
        self.assertEqual(regex_to_glob.convert(".*?\..*\.cons\.tax\.summary"), "*.*.cons.tax.summary") 
        self.assertEqual(regex_to_glob.convert(".+\.sdf$"), "*.sdf?")
    
    def test_convert_hard(self) -> None:
        self.assertEqual(regex_to_glob.convert("[^_]+_[^_]+\.fq"), "*_*.fq")
        self.assertEqual(regex_to_glob.convert("\S+\ \S+\ single 1\..*"), "* * single 1.*")
        self.assertEqual(regex_to_glob.convert(".*? index (0|1)\..*"), "* index {0,1}.*")
        self.assertEqual(regex_to_glob.convert("\S+ \S+ (single (0|1)|(forward|reverse) 0)\..*"), "* * {single {0,1},{forward,reverse} 0}.*")
 

class TestAccessoryFiles(unittest.TestCase):

    def setUp(self) -> None:
        _reset_global_settings()
        self.srcfmt = 'galaxy'
    
    def test_scripts_files_to_create(self) -> None:
        settings.translate.MODE = 'extended'
        filepath = os.path.abspath(f'{GALAXY_TESTDATA_PATH}/limma_voom_wf.ga')
        wf = ingest(filepath, self.srcfmt)
        assert(isinstance(wf, Workflow))
        tool = wf.step_nodes['limma_voom'].tool
        assert(isinstance(tool, CommandTool))
        
        # checking files_to_create entry for script
        self.assertEqual(len(tool.files_to_create()), 1)
        self.assertIn('limma_voom.R', tool.files_to_create())
    
    def test_configfiles_files_to_create(self) -> None:
        filepath = os.path.abspath(f'{GALAXY_TESTDATA_PATH}/annotate-my-ids-wf.ga')
        wf = ingest(filepath, self.srcfmt)
        assert(isinstance(wf, Workflow))
        tool = wf.step_nodes['annotatemyids'].tool
        assert(isinstance(tool, CommandTool))

        # checking files_to_create entry for configfile
        self.assertEqual(len(tool.files_to_create()), 1)
        self.assertIn('annotatemyids_script', tool.files_to_create())
    
    def test_scripts_as_params(self) -> None:
        filepath = os.path.abspath(f'{GALAXY_TESTDATA_PATH}/limma_voom_wf.ga')
        wf = ingest(filepath, self.srcfmt)
        assert(isinstance(wf, Workflow))
        tool = wf.step_nodes['limma_voom'].tool
        assert(isinstance(tool, CommandTool))

        # checking ToolInput for script
        self.assertIn('limma_voom_script', tool.inputs_map())
        tinput = [x for x in tool.inputs() if x.id() == 'limma_voom_script'][0]
        self.assertIsInstance(tinput.input_type, File)
        self.assertEqual(tinput.position, 1)
        self.assertIsNone(tinput.prefix)
    
    def test_configfiles_as_params(self) -> None:
        filepath = os.path.abspath(f'{GALAXY_TESTDATA_PATH}/annotate-my-ids-wf.ga')
        wf = ingest(filepath, self.srcfmt)
        assert(isinstance(wf, Workflow))
        tool = wf.step_nodes['annotatemyids'].tool
        assert(isinstance(tool, CommandTool))

        # checking ToolInput for configfile
        self.assertIn('annotatemyids_script', tool.inputs_map())
        tinput = [x for x in tool.inputs() if x.id() == 'annotatemyids_script'][0]
        self.assertIsInstance(tinput.input_type, File)
        self.assertEqual(tinput.position, 1)
        self.assertIsNone(tinput.prefix)
    
    def test_scripts_workflow_components(self) -> None:
        filepath = os.path.abspath(f'{GALAXY_TESTDATA_PATH}/limma_voom_wf.ga')
        wf = ingest(filepath, self.srcfmt)
        assert(isinstance(wf, Workflow))

        # checking Workflow InputNode for script
        self.assertIn('limma_voom_script', wf.input_nodes)
        
        # checking Workflow InputNode source for script when calling tool
        step = wf.step_nodes['limma_voom']
        self.assertIn('limma_voom_script', step.sources)
        source = step.sources['limma_voom_script'].source_map[0].source
        self.assertEqual(source.id(), 'limma_voom_script')
        
    def test_configfiles_workflow_components(self) -> None:
        filepath = os.path.abspath(f'{GALAXY_TESTDATA_PATH}/annotate-my-ids-wf.ga')
        wf = ingest(filepath, self.srcfmt)
        assert(isinstance(wf, Workflow))
        
        # checking Workflow InputNode for configfile
        self.assertIn('annotatemyids_script', wf.input_nodes)
        
        # checking Workflow InputNode source for configfile when calling tool
        step = wf.step_nodes['annotatemyids']
        self.assertIn('annotatemyids_script', step.sources)
        source = step.sources['annotatemyids_script'].source_map[0].source
        self.assertEqual(source.id(), 'annotatemyids_script')
    



class TestResolveDependencies(unittest.TestCase):

    def setUp(self) -> None:
        _reset_global_settings()
        settings.ingest.galaxy.GEN_IMAGES = True
        settings.ingest.galaxy.DISABLE_CONTAINER_CACHE = True
        settings.testing.TESTING_USE_DEFAULT_CONTAINER = False

    def test_coreutils_requirement(self) -> None:
        filepath = os.path.abspath(f'{GALAXY_TESTDATA_PATH}/text_processing-d698c222f354/cut.xml')
        jtool = ingest_galaxy(filepath)
        image_uri = jtool.container()
        self.assertEqual(image_uri, 'quay.io/biocontainers/coreutils:8.25--1')
    
    def test_single_requirement(self) -> None:
        wf_path = os.path.abspath(f'{GALAXY_TESTDATA_PATH}/wf_abricate.ga')
        xmltool = _load_xmltool_for_step(wf_path, 1)
        image_uri = resolve_dependencies_as_container(xmltool)
        self.assertEqual(image_uri, 'quay.io/biocontainers/abricate:1.0.1--ha8f3691_1')

    #@unittest.skipUnless(_docker_not_running(), 'docker daemon must be running to test this')
    @unittest.skip("requires docker to test")
    def test_multiple_requirements(self) -> None:
        wf_path = os.path.abspath(f'{GALAXY_TESTDATA_PATH}/hisat2_wf.ga')
        xmltool = _load_xmltool_for_step(wf_path, 2)
        image_uri = resolve_dependencies_as_container(xmltool)
        self.assertEqual(image_uri, 'quay.io/ppp-janis-translate/hisat2:2.2.1')


class TestMarkMainStatement(unittest.TestCase):

    def setUp(self) -> None:
        _reset_global_settings()

    def get_marked_main_statement(self, filepath: str, step: int) -> str:
        xmltool = _load_xmltool_for_step(filepath, step)
        text = xmltool.raw_command
        text = simplify_cmd(text, 'main_statement')
        text = mark_main_statement(text, xmltool)
        return text

    def test_wf_featurecounts(self) -> None:
        filepath = os.path.abspath(f'{GALAXY_TESTDATA_PATH}/wf_featurecounts.ga')
        text = self.get_marked_main_statement(filepath, 1)
        self.assertIn('\n__JANIS_MAIN__\n\nfeatureCounts', text)
        
    def test_wf_unicycler(self) -> None:
        filepath = os.path.abspath(f'{GALAXY_TESTDATA_PATH}/unicycler_assembly.ga')

        # fastqc
        text = self.get_marked_main_statement(filepath, 3)
        self.assertIn('\n__JANIS_MAIN__\n\nfastqc', text)
        
        # unicycler
        text = self.get_marked_main_statement(filepath, 5)
        self.assertIn('\n__JANIS_MAIN__\n\nunicycler', text)
        
        # nanoplot
        text = self.get_marked_main_statement(filepath, 6)
        self.assertIn('\n__JANIS_MAIN__\n\n$reads_temp.append("read." + str($extension))\n#end if\nNanoPlot', text)
        
        # quast
        text = self.get_marked_main_statement(filepath, 7)
        self.assertIn('echo $labels &&\n__JANIS_MAIN__\n\n#else', text)
        
        # busco
        text = self.get_marked_main_statement(filepath, 8)
        self.assertIn('\n__JANIS_MAIN__\n\n#end if\nbusco', text)

    def test_wf_rna_seq_reads_to_counts(self) -> None:
        filepath = os.path.abspath(f'{GALAXY_TESTDATA_PATH}/rna_seq_reads_to_counts.ga')

        # fastqc
        text = self.get_marked_main_statement(filepath, 2)
        self.assertIn('\n__JANIS_MAIN__\n\nfastqc', text)
        
        # cutadapt
        text = self.get_marked_main_statement(filepath, 3)
        self.assertIn('\n__JANIS_MAIN__\n\n#end if\ncutadapt', text)
        
        # hisat2
        text = self.get_marked_main_statement(filepath, 5)
        self.assertIn('\n__JANIS_MAIN__\n\n#end if\nhisat2', text)
        
        # featurecounts
        text = self.get_marked_main_statement(filepath, 6)
        self.assertIn('\n__JANIS_MAIN__\n\nfeatureCounts', text)
        
        # picard MarkDuplicates
        text = self.get_marked_main_statement(filepath, 7)
        self.assertIn('\n__JANIS_MAIN__\n\npicard\nMarkDuplicates', text)
        
        # samtools idxstats
        text = self.get_marked_main_statement(filepath, 8)
        self.assertIn('\n__JANIS_MAIN__\n\n#end if\n#end if\nsamtools idxstats', text)
        
        # rseqc geneBody_coverage
        text = self.get_marked_main_statement(filepath, 9)
        self.assertIn('\n__JANIS_MAIN__\n\ngeneBody_coverage.py', text)
        
        # rseqc infer_experiment
        text = self.get_marked_main_statement(filepath, 10)
        self.assertIn('\n__JANIS_MAIN__\ninfer_experiment.py', text)
        
        # rseqc read_distribution
        text = self.get_marked_main_statement(filepath, 11)
        self.assertIn('\n__JANIS_MAIN__\nread_distribution.py', text)
        
        # multiqc
        text = self.get_marked_main_statement(filepath, 13)
        self.assertIn('\n__JANIS_MAIN__\n\n#end for\n#end if\n#end for\nmultiqc', text)
        


class TestCommandAnnotation(unittest.TestCase):

    def setUp(self) -> None:
        _reset_global_settings()

    def test_simple_inline_bool_annotator(self) -> None:
        filepath = os.path.abspath(f'{GALAXY_TESTDATA_PATH}/cutadapt_wf.ga')
        xmltool = _load_xmltool_for_step(filepath, 2)
        command = gen_command(xmltool, annotators=['SimpleInlineBoolAnnotator'])
        expected_flags = set([
            '--no-indels', '--revcomp', '--discard-trimmed', '--discard-untrimmed',
            '--discard-cassava', '--trim-n', '--zero-cap'
        ])
        actual_flags = set(list(command.flags.keys()))
        self.assertSetEqual(actual_flags, expected_flags)
        self.assertEqual(len(command.options), 0)
        self.assertEqual(len(command.positionals), 0)
        self.assertIsNone(command.redirect)

    def test_simple_multiline_bool_annotator(self) -> None:
        filepath = os.path.abspath(f'{GALAXY_TESTDATA_PATH}/limma_voom_wf.ga')
        xmltool = _load_xmltool_for_step(filepath, 3)
        command = gen_command(xmltool, annotators=['SimpleMultilineBoolAnnotator'])
        # expected_flags = set([
        #     '-y', '-F', '-x', '-L', '-r', '-T', '-w', '-b'
        # ])
        expected_flags = set([
            '-F', '-x', '-L', '-r', '-T', '-w', '-b'
        ])
        actual_flags = set(list(command.flags.keys()))
        self.assertSetEqual(actual_flags, expected_flags)
        self.assertEqual(len(command.options), 0)
        self.assertEqual(len(command.positionals), 0)
        self.assertIsNone(command.redirect)

    def test_simple_select_annotator(self) -> None:
        filepath = os.path.abspath(f'{GALAXY_TESTDATA_PATH}/cutadapt_wf.ga')
        xmltool = _load_xmltool_for_step(filepath, 2)
        command = gen_command(xmltool, annotators=['SimpleSelectAnnotator'])
        expected_flags = set([
            '--match-read-wildcards', '--no-match-adapter-wildcards'
        ])
        expected_options = set([
            '--action', '--pair-filter'
        ])
        actual_options = set(list(command.options.keys()))
        actual_flags = set(list(command.flags.keys()))
        self.assertSetEqual(actual_flags, expected_flags)
        self.assertSetEqual(actual_options, expected_options)
        self.assertEqual(len(command.positionals), 0)
        self.assertIsNone(command.redirect)

    def test_option_param_annotator(self) -> None:
        filepath = os.path.abspath(f'{GALAXY_TESTDATA_PATH}/cutadapt_wf.ga')
        xmltool = _load_xmltool_for_step(filepath, 2)
        command = gen_command(xmltool, annotators=['OptionParamAnnotator'])
        expected_options = set([
            '--nextseq-trim',
            '-U',
            '--quality-cutoff',
            '--length',
            '--length-tag',
            '--times',
            '--overlap',
            '--rename',
            '--strip-suffix',
            '--max-expected-errors',
            '-u',
            '-Q',
            '--max-n',
            '--error-rate',
        ])
        actual_options = set(list(command.options.keys()))
        self.assertSetEqual(actual_options, expected_options)
        self.assertEqual(len(command.flags), 0)
        self.assertEqual(len(command.positionals), 0)
        self.assertIsNone(command.redirect)

    def test_local_cmdstr_annotator(self) -> None:
        filepath = os.path.abspath(f'{GALAXY_TESTDATA_PATH}/limma_voom_wf.ga')
        xmltool = _load_xmltool_for_step(filepath, 3)
        command = gen_command(xmltool, annotators=['LocalCmdstrAnnotator'])
        expected_options = set([
            '-C',
            '-t',
            '-P',
            '-z',
            '-G',
            '-a',
            '-n',
            '-c',
            '-m',
            '-l',
            '-f',
            '-d',
            '-p',
            '-s',
        ])
        actual_options = set(list(command.options.keys()))
        self.assertSetEqual(actual_options, expected_options)
        self.assertEqual(len(command.flags), 0)
        self.assertEqual(len(command.positionals), 0)
        self.assertIsNone(command.redirect)

    def test_global_cmdstr_annotator(self) -> None:
        filepath = os.path.abspath(f'{GALAXY_TESTDATA_PATH}/limma_voom_wf.ga')
        xmltool = _load_xmltool_for_step(filepath, 3)
        command = gen_command(xmltool, annotators=['GlobalCmdstrAnnotator'])
        expected_positionals = set([
            'Rscript', 'limma_voom_script'
        ])
        expected_flags = set([
            '-x',
            '-T',
            '-r',
            '-F',
            '-y',
            '-b',
            '-w',
            '-L',

        ])
        expected_options = set([
            '-i',
            '-n',
            '-P',
            '-a',
            '-m',
            '-G',
            '-t',
            '-c',
            '-f',
            '-o',
            '-C',
            '-R',
            '-s',
            '-j',
            '-l',
            '-D',
            '-z',
            '-d',
            '-p',
        ])
        actual_flags = set(list(command.flags.keys()))
        actual_options = set(list(command.options.keys()))
        actual_positionals = set([x.name for x in command.positionals.values()])
        self.assertSetEqual(actual_positionals, expected_positionals)
        self.assertSetEqual(actual_flags, expected_flags)
        self.assertSetEqual(actual_options, expected_options)
        self.assertIsNone(command.redirect)

    def test_all_goseq(self) -> None:
        filepath = os.path.abspath(f'{GALAXY_TESTDATA_PATH}/goseq_wf.ga')
        xmltool = _load_xmltool_for_step(filepath, 2)
        command = gen_command(xmltool)
        expected_positionals = set([
            'Rscript', 'goseq_script'
        ])
        expected_flags = set([])
        expected_options = set([
            '--genome',
            '--sample_vs_wallenius_plot',
            '--nobias_tab',
            '--p_adj_method',
            '--repcnt',
            '--length_file',
            '--gene_id',
            '--use_genes_without_cat',
            '--categories_genes_out_fp',
            '--dge_file',
            '--wallenius_tab',
            '--length_bias_plot',
            '--rdata',
            '--fetch_cats',
            '--make_plots',
            '--top_plot',
            '--category_file',
            '--sampling_tab',
        ])
        actual_flags = set(list(command.flags.keys()))
        actual_options = set(list(command.options.keys()))
        actual_positionals = set([x.name for x in command.positionals.values()])
        self.assertSetEqual(actual_positionals, expected_positionals)
        self.assertSetEqual(actual_flags, expected_flags)
        self.assertSetEqual(actual_options, expected_options)
        self.assertIsNone(command.redirect)
        
    def test_all_limma_voom(self) -> None:
        filepath = os.path.abspath(f'{GALAXY_TESTDATA_PATH}/limma_voom_wf.ga')
        xmltool = _load_xmltool_for_step(filepath, 3)
        command = gen_command(xmltool)
        expected_positionals = set([
            'Rscript', 'limma_voom_script'
        ])
        expected_flags = set([
            '-x',
            '-T',
            '-r',
            '-F',
            '-y',
            '-b',
            '-w',
            '-L',
            '-y',
        ])
        expected_options = set([
            '-i',
            '-n',
            '-P',
            '-a',
            '-m',
            '-G',
            '-t',
            '-c',
            '-f',
            '-o',
            '-C',
            '-R',
            '-s',
            '-j',
            '-l',
            '-D',
            '-z',
            '-d',
            '-p',
        ])
        actual_flags = set(list(command.flags.keys()))
        actual_options = set(list(command.options.keys()))
        actual_positionals = set([x.name for x in command.positionals.values()])
        self.assertSetEqual(actual_positionals, expected_positionals)
        self.assertSetEqual(actual_flags, expected_flags)
        self.assertSetEqual(actual_options, expected_options)
        self.assertIsNone(command.redirect)
    
    def test_all_cutadapt(self) -> None:
        filepath = os.path.abspath(f'{GALAXY_TESTDATA_PATH}/cutadapt_wf.ga')
        xmltool = _load_xmltool_for_step(filepath, 2)
        command = gen_command(xmltool)
        expected_positionals = set([
            'cutadapt', 'input_1', 'input_2'
        ])
        expected_flags = set([
            '--no-indels', 
            '--revcomp', 
            '--discard-trimmed', 
            '--discard-untrimmed',
            '--discard-cassava', 
            '--trim-n', 
            '--zero-cap', 
            '--match-read-wildcards', 
            '--no-match-adapter-wildcards'
        ])
        expected_options = set([
            '--action', 
            '--pair-filter',
            '-u',
            '-U',
            '-Q',
            '--error-rate',
            '--times',
            '--overlap',
            '--max-n',
            '--max-expected-errors',
            '--quality-cutoff',
            '--nextseq-trim',
            '--strip-suffix',
            '--length',
            '--length-tag',
            '--rename',
            '--minimum-length',
            '--maximum-length',
            '-j',
            '--json',
            '--info-file',
            '--rest-file',
            '--wildcard-file',
            '--too-short-output',
            '--too-long-output',
            '--too-short-paired-output',
            '--too-long-paired-output',
            '--untrimmed-paired-output',
            '--paired-output',
            '--untrimmed-output',
            '--output',
        ])
        actual_flags = set(list(command.flags.keys()))
        actual_options = set(list(command.options.keys()))
        actual_positionals = set([x.name for x in command.positionals.values()])
        self.assertSetEqual(actual_positionals, expected_positionals)
        self.assertSetEqual(actual_flags, expected_flags)
        self.assertSetEqual(actual_options, expected_options)
        self.assertIsNotNone(command.redirect)

    def test_all_hisat2(self) -> None:
        filepath = os.path.abspath(f'{GALAXY_TESTDATA_PATH}/hisat2_wf.ga')
        gx_workflow = _load_gxworkflow(filepath)
        gx_step = gx_workflow['steps']['2']
        _configure_tool_settings(gx_step)
        xmltool = load_xmltool(runtime.tool.tool_path)
        command = gen_command(xmltool, gx_step)
        expected_positionals = set([
            'hisat2'
        ])
        expected_flags = set([
            '--tmo',
            '--nofw',
            '--non-deterministic',
            '--no-softclip',
            '--dta-cufflinks',
            '--ff',
            '--no-discordant',
            '--phred64',
            '--ignore-quals',
            '--fr',
            '--no-unal',
            '--int-quals',
            '--dta',
            '--phred33',
            '--new-summary',
            '--add-chrname',
            '--omit-sec-seq',
            '--norc',
            '--no-mixed',
            '--remove-chrname',
            '--no-templatelen-adjustment',
            '--rf',
            '-f',
            '--no-spliced-alignment',
            '--solexa-quals',
        ])
        expected_options = set([
            '--pen-cansplice',
            '--mp',
            '--rdg',
            '-X',
            '--al-gz',
            '--sp',
            '--score-min',
            '--pen-noncanintronlen',
            '-I',
            '--n-ceil',
            '--trim3',
            '--trim5',
            '--seed',
            '--min-intronlen',
            '--rfg',
            '--known-splicesite-infile',
            '--un-gz',
            '--max-intronlen',
            '-2',
            '--qupto',
            '-x',
            '--skip',
            '--al-bz2',
            '-p',
            '-1',
            '-k',
            '--summary-file',
            '--np',
            '--pen-noncansplice',
            '-U',
            '--novel-splicesite-outfile',
            '--rna-strandness',
            '--un-bz2',
            '--pen-canintronlen',
            '--rg-id',
        ])
        actual_flags = set(list(command.flags.keys()))
        actual_options = set(list(command.options.keys()))
        actual_positionals = set([x.name for x in command.positionals.values()])
        self.assertSetEqual(actual_positionals, expected_positionals)
        self.assertSetEqual(actual_flags, expected_flags)
        self.assertSetEqual(actual_options, expected_options)
        self.assertIsNone(command.redirect)

    def test_all_picard_markduplicates(self) -> None:
        filepath = os.path.abspath(f'{GALAXY_TESTDATA_PATH}/wf_mark_duplicates.ga')
        gx_workflow = _load_gxworkflow(filepath)
        gx_step = gx_workflow['steps']['1']
        _configure_tool_settings(gx_step)
        xmltool = load_xmltool(runtime.tool.tool_path)
        command = gen_command(xmltool, gx_step)

        expected_positionals = set([
            'picard',
            'MarkDuplicates',
        ])

        expected_flags = set()

        expected_options = set([
            'OPTICAL_DUPLICATE_PIXEL_DISTANCE',
            'OUTPUT',
            'QUIET',
            'METRICS_FILE',
            'BARCODE_TAG',
            'REMOVE_DUPLICATES',
            'VALIDATION_STRINGENCY',
            'VERBOSITY',
            'DUPLICATE_SCORING_STRATEGY',
            'TAGGING_POLICY',
            'INPUT',
            'ASSUME_SORTED'
        ])

        actual_flags = set(list(command.flags.keys()))
        actual_options = set(list(command.options.keys()))
        actual_positionals = set([x.name for x in command.positionals.values()])

        self.assertSetEqual(actual_positionals, expected_positionals)
        self.assertSetEqual(actual_flags, expected_flags)
        self.assertSetEqual(actual_options, expected_options)
        self.assertIsNone(command.redirect)
    
    def test_all_samtools_idxstats(self) -> None:
        filepath = os.path.abspath(f'{GALAXY_TESTDATA_PATH}/wf_samtools_idxstats.ga')
        gx_workflow = _load_gxworkflow(filepath)
        gx_step = gx_workflow['steps']['1']
        _configure_tool_settings(gx_step)
        xmltool = load_xmltool(runtime.tool.tool_path)
        command = gen_command(xmltool, gx_step)

        expected_positionals = set([
            'samtools',
            'idxstats',
            'input'
        ])
        expected_flags = set()
        expected_options = set([
            '-@'
        ])

        actual_flags = set(list(command.flags.keys()))
        actual_options = set(list(command.options.keys()))
        actual_positionals = set([x.name for x in command.positionals.values()])

        self.assertSetEqual(actual_positionals, expected_positionals)
        self.assertSetEqual(actual_flags, expected_flags)
        self.assertSetEqual(actual_options, expected_options)
        self.assertIsNotNone(command.redirect)



class TestLoadToolState(unittest.TestCase):
    """
    tests ability to map or generate janis datatypes / selectors etc
    from internal tool objects.
    """
    def setUp(self) -> None:
        _reset_global_settings()
        

    def test_default_filters(self) -> None:
        wf_path = os.path.abspath(f'{GALAXY_TESTDATA_PATH}/cutadapt_wf.ga')
        gx_workflow = _load_gxworkflow(wf_path)
        gx_step = gx_workflow['steps']['2']
        tool_state = _load_tool_state(gx_step)
        self.assertEqual(tool_state['output_selector'][0], 'report')
        self.assertEqual(tool_state['adapter_options']['action'], 'trim')
        self.assertEqual(tool_state['library']['r1']['cut'], '0')
    
    def test_null_varname_filter(self) -> None:
        wf_path = os.path.abspath(f'{GALAXY_TESTDATA_PATH}/cutadapt_wf.ga')
        gx_workflow = _load_gxworkflow(wf_path)
        gx_step = gx_workflow['steps']['2']
        tool_state = _load_tool_state(gx_step, additional_filters=['ReplaceNullWithVarname'])
        self.assertEqual(tool_state['filter_options']['minimum_length'], '$filter_options.minimum_length')
        self.assertEqual(tool_state['filter_options']['maximum_length'], '$filter_options.maximum_length')
        self.assertEqual(tool_state['filter_options']['max_n'], '$filter_options.max_n')
        self.assertEqual(tool_state['filter_options']['max_expected_errors'], '$filter_options.max_expected_errors')
        self.assertEqual(tool_state['adapter_options']['action'], 'trim')
        self.assertEqual(tool_state['library']['r1']['cut'], '0')

    def test_flat_filter(self) -> None:
        wf_path = os.path.abspath(f'{GALAXY_TESTDATA_PATH}/cutadapt_wf.ga')
        gx_workflow = _load_gxworkflow(wf_path)
        gx_step = gx_workflow['steps']['2']
        tool_state = _load_tool_state(gx_step, additional_filters=['Flatten'])
        self.assertEqual(tool_state['output_selector'][0], 'report')
        self.assertEqual(tool_state['adapter_options.action'], 'trim')
        self.assertEqual(tool_state['filter_options.minimum_length'], None)
        self.assertEqual(tool_state['library.r1.cut'], '0')




class TestJanisGeneralMapping(unittest.TestCase):
    """
    tests ability to map or generate janis datatypes / selectors etc
    from internal tool objects.
    """
    def setUp(self) -> None:
        _reset_global_settings()
        self.tool = MOCK_TOOL_ABRICATE

    def test_to_janis_datatype(self) -> None:
        """
        checks that JanisDatatype objects can be mapped to
        janis DataType objects 
        """
        jtype1 = to_janis_datatype(self.tool.inputs[0])  
        jtype2 = to_janis_datatype(self.tool.inputs[1])
        jtype3 = to_janis_datatype(self.tool.inputs[2])
        jtype4 = to_janis_datatype(self.tool.outputs[0])
        # check janis objects are correct
        self.assertIsInstance(jtype1, File)
        self.assertIsInstance(jtype2, Boolean)
        self.assertIsInstance(jtype3, Array)
        self.assertIsInstance(jtype4, Stdout)
        # check attributes are correct
        self.assertEquals(jtype2.optional, True)
        self.assertIsInstance(jtype3.subtype(), Float)
        self.assertIsInstance(jtype4.subtype, File)

    def test_to_janis_selector(self) -> None:
        """
        checks that OutputComponent objects can generate 
        janis Selector objects
        """
        jsel1 = to_janis_selector(self.tool.outputs[0])
        jsel2 = to_janis_selector(self.tool.outputs[1])
        jsel3 = to_janis_selector(self.tool.outputs[2])
        # check janis objects are correct
        self.assertIsInstance(jsel2, WildcardSelector)
        self.assertIsInstance(jsel3, InputSelector)
        self.assertIsNone(jsel1)
        # check attributes are correct
        self.assertEquals(jsel2.wildcard, 'report.txt')
        self.assertEquals(jsel3.input_to_select, 'file_input')



class TestJanisToolMapping(unittest.TestCase):
    """
    tests ability to map or generate janis tool objects 
    from internal tool objects.
    MOCK_TOOL is used to test mapping functions.
    """
    def setUp(self) -> None:
        _reset_global_settings()
        self.tool = MOCK_TOOL_ABRICATE
    
    def test_to_janis_tool(self) -> None:
        """
        checks that Tool objects can be mapped to
        janis CommandToolBuilder objects 
        """
        jtool = to_janis_tool(self.tool)
        self.assertIsInstance(jtool, CommandToolBuilder)

    def test_to_janis_tool_input(self) -> None:
        """
        checks that InputComponent objects can be mapped to
        janis ToolInput objects 
        """
        jinp1 = to_janis_tool_input(self.tool.inputs[0])  
        jinp2 = to_janis_tool_input(self.tool.inputs[1])
        jinp3 = to_janis_tool_input(self.tool.inputs[2])
        jinp4 = to_janis_tool_input(self.tool.inputs[3])
        # check janis objects are correct
        self.assertIsInstance(jinp1, ToolInput)
        self.assertIsInstance(jinp2, ToolInput)
        self.assertIsInstance(jinp3, ToolInput)
        self.assertIsInstance(jinp4, ToolInput)
        # check attributes are correct
        self.assertEquals(jinp1.tag, 'file_input')
        self.assertEquals(jinp1.prefix, None)
        self.assertEquals(jinp1.separate_value_from_prefix, True)
        self.assertIsInstance(jinp1.input_type, File)
        
        self.assertEquals(jinp2.tag, 'noheader')
        self.assertEquals(jinp2.prefix, '--noheader')
        self.assertIsInstance(jinp2.input_type, Boolean)
        self.assertEquals(jinp2.input_type.optional, True)

        self.assertEquals(jinp3.tag, 'adv_min_dna_id')
        self.assertEquals(jinp3.prefix, '--minid=')
        self.assertEquals(jinp3.separate_value_from_prefix, False)
        self.assertIsInstance(jinp3.input_type, Array)
        
        self.assertEquals(jinp4.tag, 'adv_db')
        self.assertEquals(jinp4.prefix, '--db=')
        self.assertEquals(jinp4.separate_value_from_prefix, False)
        self.assertEquals(jinp4.default, 'resfinder')
        self.assertIsInstance(jinp4.input_type, String)

    def test_to_janis_tool_output(self) -> None:
        """
        checks that OutputComponent objects can be mapped to
        janis ToolOutput objects 
        """
        jout1 = to_janis_tool_output(self.tool.outputs[0])
        jout2 = to_janis_tool_output(self.tool.outputs[1])
        jout3 = to_janis_tool_output(self.tool.outputs[2])
        # check janis objects are correct
        self.assertIsInstance(jout1, ToolOutput)
        self.assertIsInstance(jout2, ToolOutput)
        self.assertIsInstance(jout3, ToolOutput)
        # check attributes are correct
        self.assertIsInstance(jout1.output_type, Stdout)
        self.assertEquals(jout1.tag, 'out_report1')
        self.assertIsNone(jout1.selector)
        
        self.assertIsInstance(jout2.output_type, File)
        self.assertEquals(jout2.tag, 'out_report2')
        self.assertIsInstance(jout2.selector, WildcardSelector)
        self.assertEquals(jout2.selector.wildcard, 'report.txt')
        
        self.assertIsInstance(jout3.output_type, File)
        self.assertEquals(jout3.tag, 'out_file_input')
        self.assertIsInstance(jout3.selector, InputSelector)
        self.assertEquals(jout3.selector.input_to_select, 'file_input')

    def test_to_janis_metadata(self) -> None:
        """
        checks that ToolMetadata objects can be mapped to
        janis ToolXMLMetadata objects 
        """
        # MOCK_TOOL
        jmeta = to_janis_metadata(self.tool.metadata)
        self.assertIsInstance(jmeta, ToolMetadata)
        self.assertIsNotNone(jmeta.contributors)
        self.assertIsNotNone(jmeta.dateCreated)
        self.assertIsNotNone(jmeta.dateUpdated)
        self.assertIsNotNone(jmeta.citation)
        self.assertIsNotNone(jmeta.documentation)
        self.assertIsNotNone(jmeta.short_documentation)
        self.assertIsNotNone(jmeta.version)


class TestJanisWorkflowMapping(unittest.TestCase):
    """
    tests ability to map or generate janis workflow objects 
    from internal galaxy ingest model
    """
    def setUp(self) -> None:
        _reset_global_settings()
        self.internal = MOCK_WORKFLOW
        self.jworkflow = to_janis_workflow(self.internal)
        self.jinputs = to_janis_inputs_dict(self.internal)

    def test_to_janis_workflow(self) -> None:
        self.assertIsInstance(self.jworkflow, WorkflowBuilder)
    
    def test_to_janis_inputs_dict(self) -> None:
        # single input, no value
        self.assertEquals(self.jinputs['in_fasta'], None)
        # single input, provided value
        self.internal.inputs[0].value = 'path/to/file.fasta'
        jinputs = to_janis_inputs_dict(self.internal)
        self.assertEquals(jinputs['in_fasta'], 'path/to/file.fasta')

    def test_janis_metadata(self) -> None:
        self.assertIsInstance(self.jworkflow.metadata, WorkflowMetadata)

    def test_janis_inputs(self) -> None:
        target_inputs = {
            'in_fasta', 
            'abricate_adv_db',
            'abricate_adv_min_dna_id',
            'abricate_noheader',
        }
        actual_inputs = set(self.jworkflow.input_nodes.keys())
        self.assertEquals(target_inputs, actual_inputs)

        jinp = self.jworkflow.input_nodes['in_fasta']
        self.assertIsInstance(jinp, InputNode)
        self.assertIsInstance(jinp.datatype, File)

    def test_janis_outputs(self) -> None:
        target_outputs = {'abricate_out_report1'}
        actual_outputs = set(self.jworkflow.output_nodes.keys())
        self.assertEquals(target_outputs, actual_outputs)

        jout = self.jworkflow.output_nodes['abricate_out_report1']
        self.assertIsInstance(jout, OutputNode)
        self.assertIsInstance(jout.datatype, Stdout)
        self.assertIsInstance(jout.datatype.subtype, File)
        self.assertEquals(jout.doc.doc, 'report file')

    def test_janis_steps(self) -> None:
        target_steps = {'abricate'}
        actual_steps = set(self.jworkflow.step_nodes.keys())
        self.assertEquals(target_steps, actual_steps)

        # basic object checks
        jstep = self.jworkflow.step_nodes['abricate']
        self.assertIsInstance(jstep, StepNode)
        self.assertIsInstance(jstep.tool, CommandTool)

        # sources (tool input values)
        self.assertIn('file_input', jstep.sources)
        self.assertIn('noheader', jstep.sources)
        self.assertIn('adv_min_dna_id', jstep.sources)
        self.assertIn('adv_db', jstep.sources)
        
        # scatter 
        self.assertEquals(jstep.scatter.fields, ['adv_min_dna_id'])
        self.assertEquals(jstep.scatter.method, ScatterMethod.dot)



class TestDatatypeInference(unittest.TestCase):
    """
    tests the datatype which is assigned to an entity
    """
    def setUp(self) -> None:
        _reset_global_settings()

    def test_positional(self) -> None:
        dtype = datatypes.get(MOCK_POSITIONAL1)
        self.assertEquals(dtype.classname, 'Fastq')
     
    def test_flag(self) -> None:
        self.assertEquals(datatypes.get(MOCK_FLAG1), bool_t)
    
    def test_option(self) -> None:
        self.assertEquals(datatypes.get(MOCK_OPTION2), string_t)
    
    def test_outputs(self) -> None:
        dtype = datatypes.get(MOCK_REDIRECT_OUTPUT)
        self.assertEquals(dtype.classname, 'TextFile')
    
    def test_workflow_input(self) -> None:
        self.assertEquals(datatypes.get(MOCK_WORKFLOW_INPUT1), file_t)


class TestFromGalaxy(unittest.TestCase):

    def setUp(self) -> None:
        _reset_global_settings()

    def test_ingest_abricate_tool(self) -> None:
        filepath = os.path.abspath(f'{GALAXY_TESTDATA_PATH}/abricate/abricate.xml')
        jtool = ingest_galaxy(filepath)
        assert(isinstance(jtool, CommandTool))
        
        self.assertEquals(len(jtool.inputs()), 5)
        self.assertEquals(len(jtool.outputs()), 1)
        self.assertEquals(jtool.base_command(), ['abricate'])

    def test_ingest_cutadapt_wf(self) -> None:
        filepath = os.path.abspath(f'{GALAXY_TESTDATA_PATH}/cutadapt_wf.ga')
        jworkflow = ingest_galaxy(filepath)
        assert(isinstance(jworkflow, WorkflowBuilder))

        self.assertEquals(len(jworkflow.input_nodes), 42)
        self.assertIn('in_forward', jworkflow.input_nodes)
        self.assertIn('in_reverse', jworkflow.input_nodes)
        self.assertEquals(len(jworkflow.step_nodes), 1)
        self.assertIn('cutadapt', jworkflow.step_nodes)
        self.assertEquals(len(jworkflow.output_nodes), 3)
        self.assertIn('cutadapt_out12', jworkflow.output_nodes)
        self.assertIn('cutadapt_out22', jworkflow.output_nodes)
        self.assertIn('cutadapt_out_report', jworkflow.output_nodes)
    
    def test_ingest_unicycler_assembly(self) -> None:
        filepath = os.path.abspath(f'{GALAXY_TESTDATA_PATH}/unicycler_assembly.ga')
        jworkflow = ingest_galaxy(filepath)
        assert(isinstance(jworkflow, WorkflowBuilder))

        self.assertEquals(len(jworkflow.step_nodes), 6)
        self.assertEquals(len(jworkflow.output_nodes), 7)
        self.assertIn('in_short_R1', jworkflow.input_nodes)
        self.assertIn('in_short_R2', jworkflow.input_nodes)
        self.assertIn('in_long', jworkflow.input_nodes)
    
    def test_translate_cutadapt_wf_nextflow(self) -> None:
        srcfmt = 'galaxy'
        destfmt = 'nextflow'
        filepath = os.path.abspath(f'{GALAXY_TESTDATA_PATH}/cutadapt_wf.ga')
        _run(filepath, srcfmt, destfmt)
    
    def test_translate_unicycler_wf_nextflow(self) -> None:
        srcfmt = 'galaxy'
        destfmt = 'nextflow'
        filepath = os.path.abspath(f'{GALAXY_TESTDATA_PATH}/unicycler_assembly.ga')
        _run(filepath, srcfmt, destfmt)

