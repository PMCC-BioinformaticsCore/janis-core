

import unittest
from typing import Any

from janis_core import settings
from janis_core.ingestion import ingest
from janis_core.translations import translate
from janis_core.tests.testtools import FileOutputPythonTestTool
from janis_core.tests.testtools import GridssTestTool
from janis_core.tests.testworkflows import PruneFlatTW
from janis_core.tests.testworkflows import PruneNestedTW
from janis_core.tests.testworkflows import AssemblyTestWF
from janis_core.redefinitions.tools import Cat
from janis_core.redefinitions.tools import GenerateVardictHeaderLines
from janis_core.redefinitions.workflows import BwaAligner
from janis_core.redefinitions.workflows import WGSGermlineMultiCallers

from janis_core import CommandToolBuilder
from janis_core import WorkflowBuilder
from janis_core import CodeTool

from janis_core.translations import nextflow
from janis_core.translations.common import to_builders
from janis_core.translations.common import prune_workflow
from janis_core import settings

import os 
import regex as re
import yaml

CWL_TESTDATA_PATH = os.path.join(os.getcwd(), 'janis_core/tests/data/cwl')
GALAXY_TESTDATA_PATH = os.path.join(os.getcwd(), 'janis_core/tests/data/galaxy')
JANIS_TESTDATA_PATH = os.path.join(os.getcwd(), 'janis_core/tests/data/janis')
WDL_TESTDATA_PATH = os.path.join(os.getcwd(), 'janis_core/tests/data/wdl')


# ------- HELPER FUNCS ------- #

def _run(filepath: str, srcfmt: str, destfmt: str) -> Any:
    internal = ingest(filepath, srcfmt)
    return translate(internal, destfmt, export_path='./translated')

def _is_nf_process(filecontents: str) -> bool:
    pattern = r'process.*?\{'
    if re.findall(pattern, filecontents):
        return True
    return False

def _is_cwl_clt(filecontents: str) -> bool:
    if 'class: CommandLineTool' in filecontents:
        return True
    return False

def _is_wdl_task(filecontents: str) -> bool:
    pattern = r'task.*?\{'
    if re.findall(pattern, filecontents):
        return True
    return False

def _get_cwl_clt_inputs(clt_text: str) -> list[str]:
    spec = yaml.safe_load(clt_text)
    return spec['inputs']

def _get_cwl_clt_args(clt_text: str) -> list[str]:
    spec = yaml.safe_load(clt_text)
    return spec['arguments']

def _get_wdl_task_command_lines(task_text: str) -> list[str]:
    """Returns the lines of the process script"""
    out: list[str] = []
    lines = task_text.split('\n')
    within_script: bool = False
    
    for i in range(len(lines)):
        if lines[i].strip() == 'command <<<':
            within_script = True
            continue
        if lines[i].strip() == '>>>' and within_script:
            within_script = False
            continue
        if within_script:
            out.append(lines[i])
    
    return out

def _get_nf_process_input_lines(process_text: str) -> list[str]:
    """Returns the lines of the process script"""
    out: list[str] = []
    lines = process_text.split('\n')
    within_inputs: bool = False
    
    for i in range(len(lines)):
        if lines[i].strip() == 'input:':
            within_inputs = True
            continue
        if lines[i].strip() == 'output:' and within_inputs:
            within_inputs = False
            continue
        if within_inputs and lines[i].strip() != '':
            out.append(lines[i])
    
    return out

def _get_nf_process_script_lines(process_text: str) -> list[str]:
    """Returns the lines of the process script"""
    out: list[str] = []
    lines = process_text.split('\n')
    within_script: bool = False
    
    for i in range(len(lines)):
        if lines[i].strip() == '"""' and not within_script:
            within_script = True
            continue
        if lines[i].strip() == '"""' and within_script:
            within_script = False
            continue
        if within_script:
            out.append(lines[i])
    
    return out

def _reset_global_settings() -> None:
    nextflow.task_inputs.clear()
    nextflow.params.clear()
    settings.ingest.SAFE_MODE = False
    settings.ingest.galaxy.GEN_IMAGES = False
    settings.ingest.galaxy.DISABLE_CONTAINER_CACHE = False
    settings.ingest.cwl.INGEST_JAVASCRIPT_EXPRESSIONS = True
    settings.ingest.cwl.REQUIRE_CWL_VERSION = False
    settings.datatypes.ALLOW_UNPARSEABLE_DATATYPES = True
    settings.graph.ALLOW_UNKNOWN_SOURCE = True
    settings.graph.ALLOW_UNKNOWN_SCATTER_FIELDS = True
    settings.graph.ALLOW_INCORRECT_NUMBER_OF_SOURCES = True
    settings.graph.ALLOW_NON_ARRAY_SCATTER_INPUT = True
    settings.graph.ALLOW_INCOMPATIBLE_TYPES = True
    settings.validation.STRICT_IDENTIFIERS = False
    settings.validation.VALIDATE_STRINGFORMATTERS = False
    settings.translate.MODE = 'extended'
    settings.translate.ALLOW_EMPTY_CONTAINER = True
    settings.translate.nextflow.ENTITY = 'workflow'

# ensuring janis pipelines for janis translate tests
def janis_pipelines_installed() -> bool:
    try:
        import janis_pipelines
        import janis_bioinformatics
        return True
    except:
        return False
    


# ------- TRANSLATION ENDPOINTS ------- #

class TestTranslationEndpoints(unittest.TestCase):
    """
    this class tests all translation endpoints. 
    exists so we know that the endpoints are working.
    translation is always to nextflow as we're just testing endpoints, 
    and want the nextflow translation to be the best translation unit so pick this one.
    """

    def setUp(self) -> None:
        _reset_global_settings()

    # CWL INGEST -> TRANSLATE

    def test_from_cwl_to_nextflow_commandtool(self) -> None:
        filepath = f'{CWL_TESTDATA_PATH}/workflows/analysis-workflows/tools/gatk_haplotype_caller.cwl'
        _run(filepath, 'cwl', 'nextflow')

    def test_from_cwl_to_nextflow_etool(self) -> None:
        filepath = f'{CWL_TESTDATA_PATH}/workflows/ebi-metagenomics/utils/return_directory.cwl'
        _run(filepath, 'cwl', 'nextflow')

    def test_from_cwl_to_nextflow_workflow(self) -> None:
        filepath = f'{CWL_TESTDATA_PATH}/workflows/analysis-workflows/subworkflows/align_sort_markdup.cwl'
        _run(filepath, 'cwl', 'nextflow')
    
    # GALAXY INGEST -> TRANSLATE

    def test_from_galaxy_to_nextflow_tool(self) -> None:
        filepath = f'{GALAXY_TESTDATA_PATH}/abricate/abricate.xml'
        _run(filepath, 'galaxy', 'nextflow')
    
    def test_from_galaxy_to_nextflow_tool_toolshed(self) -> None:
        filepath = f'toolshed.g2.bx.psu.edu/repos/devteam/samtools_flagstat/samtools_flagstat/2.0.4'
        _run(filepath, 'galaxy', 'nextflow')

    def test_from_galaxy_to_nextflow_workflow(self) -> None:
        filepath = f'{GALAXY_TESTDATA_PATH}/cutadapt_wf.ga'
        _run(filepath, 'galaxy', 'nextflow')
    
    # WDL INGEST -> TRANSLATE

    @unittest.skip('wdl ingest needs work')
    def test_from_wdl_to_nextflow_tool(self) -> None:
        filepath = f'{WDL_TESTDATA_PATH}/bwa.wdl'
        _run(filepath, 'wdl', 'nextflow')

    @unittest.skip('wdl ingest needs work')
    def test_from_wdl_to_nextflow_workflow(self) -> None:
        filepath = f'{WDL_TESTDATA_PATH}/Reads2Map/pipelines/PreprocessingReads/PreprocessingReads.wdl'
        _run(filepath, 'wdl', 'nextflow')
    
    # JANIS.TRANSLATE()

    def test_commandtool_translate_nextflow(self) -> None:
        GridssTestTool().translate('nextflow', export_path='./translated')

    def test_codetool_translate_nextflow(self) -> None:
        FileOutputPythonTestTool().translate('nextflow', export_path='./translated')

    def test_workflow_translate_nextflow(self) -> None:
        AssemblyTestWF().translate('nextflow', export_path='./translated')

    def test_str_tool(self):
        BwaAligner().translate("janis")

    def test_str_python_tool(self):
        GenerateVardictHeaderLines().translate("janis")

    def test_command_tool(self):
        Cat().translate("janis")

    def test_str_big_workflow(self):
        WGSGermlineMultiCallers().translate("janis")



### ----- JANIS -> NEXTFLOW ----- ###

class TestJanisToNextflow(unittest.TestCase):
    
    def setUp(self) -> None:
        self.src = 'janis'
        self.dest = 'nextflow'
        _reset_global_settings()

    @unittest.skipUnless(janis_pipelines_installed(), "janis_pipelines not installed")
    def test_wf_bwa_aligner(self) -> None:
        # settings.translate.MODE = 'extended'
        # settings.translate.MODE = 'regular'
        # settings.translate.MODE = 'skeleton'
        from janis_bioinformatics.tools.common import BwaAligner
        wf = BwaAligner()
        maintask, config, subtasks = translate(wf, self.dest, export_path='./translated')

    @unittest.skipUnless(janis_pipelines_installed(), "janis_pipelines not installed")
    def test_wf_bwa_alignment(self) -> None:
        # settings.translate.MODE = 'extended'
        settings.translate.MODE = 'regular'
        # settings.translate.MODE = 'skeleton'
        from janis_pipelines.alignment.alignment import BwaAlignment
        wf = BwaAlignment()
        maintask, config, subtasks = translate(wf, self.dest, export_path='./translated')

    @unittest.skipUnless(janis_pipelines_installed(), "janis_pipelines not installed")
    def test_wf_bwa_alignment_and_qc(self) -> None:
        from janis_pipelines import BwaAlignmentAndQC
        wf = BwaAlignmentAndQC()
        maintask, config, subtasks = translate(wf, self.dest, export_path='./translated')

    @unittest.skipUnless(janis_pipelines_installed(), "janis_pipelines not installed")
    def test_wf_wgs_germline_gatk(self) -> None:
        from janis_pipelines import WGSGermlineGATK
        wf = WGSGermlineGATK()
        maintask, config, subtasks = translate(wf, self.dest, export_path='./translated')

    @unittest.skipUnless(janis_pipelines_installed(), "janis_pipelines not installed")
    def test_wf_wgs_germline_gatk_variants_only(self) -> None:
        from janis_pipelines import WGSGermlineGATKVariantsOnly
        wf = WGSGermlineGATKVariantsOnly()
        maintask, config, subtasks = translate(wf, self.dest, export_path='./translated')

    @unittest.skipUnless(janis_pipelines_installed(), "janis_pipelines not installed")
    def test_wf_wgs_germline_multi_callers(self) -> None:
        from janis_pipelines import WGSGermlineMultiCallers
        wf = WGSGermlineMultiCallers()
        maintask, config, subtasks = translate(wf, self.dest, export_path='./translated')

    @unittest.skipUnless(janis_pipelines_installed(), "janis_pipelines not installed")
    def test_wf_wgs_germline_multi_callers_variants_only(self) -> None:
        from janis_pipelines import WGSGermlineMultiCallersVariantsOnly
        wf = WGSGermlineMultiCallersVariantsOnly()
        maintask, config, subtasks = translate(wf, self.dest, export_path='./translated')

    @unittest.skipUnless(janis_pipelines_installed(), "janis_pipelines not installed")
    def test_wf_wgs_somatic_gatk(self) -> None:
        from janis_pipelines import WGSSomaticGATK
        wf = WGSSomaticGATK()
        maintask, config, subtasks = translate(wf, self.dest, export_path='./translated')

    @unittest.skipUnless(janis_pipelines_installed(), "janis_pipelines not installed")
    def test_wf_wgs_somatic_gatk_variants_only(self) -> None:
        from janis_pipelines import WGSSomaticGATKVariantsOnly
        wf = WGSSomaticGATKVariantsOnly()
        maintask, config, subtasks = translate(wf, self.dest, export_path='./translated')

    @unittest.skipUnless(janis_pipelines_installed(), "janis_pipelines not installed")
    def test_wf_wgs_somatic_multi_callers(self) -> None:
        from janis_pipelines import WGSSomaticMultiCallers
        wf = WGSSomaticMultiCallers()
        maintask, config, subtasks = translate(wf, self.dest, export_path='./translated')

    @unittest.skipUnless(janis_pipelines_installed(), "janis_pipelines not installed")
    def test_wf_wgs_somatic_multi_callers_variants_only(self) -> None:
        from janis_pipelines import WGSSomaticMultiCallersVariantsOnly
        wf = WGSSomaticMultiCallersVariantsOnly()
        maintask, config, subtasks = translate(wf, self.dest, export_path='./translated')


### ----- JANIS -> WDL ----- ###

class TestJanisToWDL(unittest.TestCase):
    
    def setUp(self) -> None:
        self.src = 'janis'
        self.dest = 'wdl'
        _reset_global_settings()

    @unittest.skipUnless(janis_pipelines_installed(), "janis_pipelines not installed")
    def test_wf_bwa_aligner(self) -> None:
        from janis_bioinformatics.tools.common import BwaAligner
        wf = BwaAligner()
        maintask, config, subtasks = translate(wf, self.dest, export_path='./translated')

    @unittest.skipUnless(janis_pipelines_installed(), "janis_pipelines not installed")
    def test_wf_bwa_alignment(self) -> None:
        from janis_pipelines.alignment.alignment import BwaAlignment
        wf = BwaAlignment()
        maintask, config, subtasks = translate(wf, self.dest, export_path='./translated')

    @unittest.skipUnless(janis_pipelines_installed(), "janis_pipelines not installed")
    def test_wf_bwa_alignment_and_qc(self) -> None:
        from janis_pipelines import BwaAlignmentAndQC
        wf = BwaAlignmentAndQC()
        maintask, config, subtasks = translate(wf, self.dest, export_path='./translated')

    @unittest.skipUnless(janis_pipelines_installed(), "janis_pipelines not installed")
    def test_wf_wgs_germline_gatk(self) -> None:
        from janis_pipelines import WGSGermlineGATK
        wf = WGSGermlineGATK()
        maintask, config, subtasks = translate(wf, self.dest, export_path='./translated')

    @unittest.skipUnless(janis_pipelines_installed(), "janis_pipelines not installed")
    def test_wf_wgs_germline_gatk_variants_only(self) -> None:
        from janis_pipelines import WGSGermlineGATKVariantsOnly
        wf = WGSGermlineGATKVariantsOnly()
        maintask, config, subtasks = translate(wf, self.dest, export_path='./translated')

    @unittest.skipUnless(janis_pipelines_installed(), "janis_pipelines not installed")
    def test_wf_wgs_germline_multi_callers(self) -> None:
        from janis_pipelines import WGSGermlineMultiCallers
        wf = WGSGermlineMultiCallers()
        maintask, config, subtasks = translate(wf, self.dest, export_path='./translated')

    @unittest.skipUnless(janis_pipelines_installed(), "janis_pipelines not installed")
    def test_wf_wgs_germline_multi_callers_variants_only(self) -> None:
        from janis_pipelines import WGSGermlineMultiCallersVariantsOnly
        wf = WGSGermlineMultiCallersVariantsOnly()
        maintask, config, subtasks = translate(wf, self.dest, export_path='./translated')

    @unittest.skipUnless(janis_pipelines_installed(), "janis_pipelines not installed")
    def test_wf_wgs_somatic_gatk(self) -> None:
        from janis_pipelines import WGSSomaticGATK
        wf = WGSSomaticGATK()
        maintask, config, subtasks = translate(wf, self.dest, export_path='./translated')

    @unittest.skipUnless(janis_pipelines_installed(), "janis_pipelines not installed")
    def test_wf_wgs_somatic_gatk_variants_only(self) -> None:
        from janis_pipelines import WGSSomaticGATKVariantsOnly
        wf = WGSSomaticGATKVariantsOnly()
        maintask, config, subtasks = translate(wf, self.dest, export_path='./translated')

    @unittest.skipUnless(janis_pipelines_installed(), "janis_pipelines not installed")
    def test_wf_wgs_somatic_multi_callers(self) -> None:
        from janis_pipelines import WGSSomaticMultiCallers
        wf = WGSSomaticMultiCallers()
        maintask, config, subtasks = translate(wf, self.dest, export_path='./translated')

    @unittest.skipUnless(janis_pipelines_installed(), "janis_pipelines not installed")
    def test_wf_wgs_somatic_multi_callers_variants_only(self) -> None:
        from janis_pipelines import WGSSomaticMultiCallersVariantsOnly
        wf = WGSSomaticMultiCallersVariantsOnly()
        maintask, config, subtasks = translate(wf, self.dest, export_path='./translated')


class TestGalaxyToNextflow(unittest.TestCase):
    
    def setUp(self) -> None:
        self.src = 'galaxy'
        self.dest = 'nextflow'
        _reset_global_settings()
        settings.translate.MODE = 'skeleton'

    # tools
    def test_tool_fastqc(self):
        filepath = os.path.abspath(f'{GALAXY_TESTDATA_PATH}/fastqc-5ec9f6bceaee/rgFastQC.xml')
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    def test_tool_nanoplot(self):
        filepath = os.path.abspath(f'{GALAXY_TESTDATA_PATH}/nanoplot-edbb6c5028f5/nanoplot.xml')
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    def test_tool_hisat2(self):
        filepath = os.path.abspath(f'{GALAXY_TESTDATA_PATH}/hisat2-6c19daec423d/hisat2.xml')
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)

    def test_tool_limma_voom_toolshed(self):
        uri = 'toolshed.g2.bx.psu.edu/repos/iuc/limma_voom/limma_voom/3.50.1+galaxy0'
        mainstr = _run(uri, self.src, self.dest)
        print(mainstr)
    
    def test_tool_samtools_flagstat_toolshed(self):
        settings.translate.MODE = 'extended'
        uri = 'toolshed.g2.bx.psu.edu/repos/devteam/samtools_flagstat/samtools_flagstat/2.0.4'
        mainstr = _run(uri, self.src, self.dest)
        print(mainstr)

    # workflows
    def test_wf_cutadapt(self):
        settings.translate.MODE = 'skeleton'
        filepath = os.path.abspath(f'{GALAXY_TESTDATA_PATH}/cutadapt_wf.ga')
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)

    def test_wf_nanoplot(self):
        filepath = os.path.abspath(f'{GALAXY_TESTDATA_PATH}/nanoplot_wf.ga')
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    def test_wf_fastqc(self):
        filepath = os.path.abspath(f'{GALAXY_TESTDATA_PATH}/fastqc_wf.ga')
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)

    def test_wf_rna_seq_reads_to_counts(self):
        settings.translate.MODE = 'skeleton'
        # settings.translate.MODE = 'regular'
        filepath = f'{GALAXY_TESTDATA_PATH}/rna_seq_reads_to_counts.ga'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)

    def test_wf_unicycler_assembly(self):
        filepath = f'{GALAXY_TESTDATA_PATH}/unicycler_assembly.ga'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    # def test_rna_seq_counts_to_genes(self):
    #     filepath = f'{GALAXY_TESTDATA_PATH}/rna_seq_counts_to_genes.ga'
    #     mainstr = _run(filepath, self.src, self.dest)
    #     print(mainstr)
    
    # def test_rna_seq_genes_to_pathways(self):
    #     filepath = f'{GALAXY_TESTDATA_PATH}/rna_seq_genes_to_pathways.ga'
    #     mainstr = _run(filepath, self.src, self.dest)
    #     print(mainstr)




# ---- PREPROCESSING: PRUNE ------------------------------

class TestPreprocessingPrune(unittest.TestCase):
    
    def setUp(self) -> None:
        _reset_global_settings()
        self.flat_wf = PruneFlatTW()
        self.flat_wf = to_builders(self.flat_wf)
        prune_workflow(self.flat_wf)
        self.nested_wf = PruneNestedTW()
        self.nested_wf = to_builders(self.nested_wf)
        prune_workflow(self.nested_wf)

    ### workflow inputs ###
    
    def test_flat_wf_inputs(self) -> None:
        actual_inputs = set(list(self.flat_wf.input_nodes.keys()))
        expected_inputs = set([
            'inFile1',
            'inFile2',
            'inFile3',
            
            'inStr1',
            'inStr2',
            'inStr3',
            
            'inFileOpt1',
            'inFileOpt2',
            
            'inStrOpt1',
            'inStrOpt2',

        ])
        self.assertSetEqual(actual_inputs, expected_inputs)
    
    def test_nested_wf_inputs(self) -> None:
        actual_inputs = set(list(self.nested_wf.input_nodes.keys()))
        expected_inputs = set([
            'inFile1',
            'inFile2',
            'inFile3',
            
            'inStr1',
            'inStr2',
            'inStr3',
            
            'inFileOpt1',
            'inFileOpt2',
            'inFileOpt3',
            
            'inStrOpt1',
            'inStrOpt2',
            'inStrOpt3',
        ])
        self.assertSetEqual(actual_inputs, expected_inputs)
        
    ### migrating statics ###

    def test_migrate_single_statics(self) -> None:
        tool = self.flat_wf.step_nodes['stp4'].tool
        
        # checking num inputs
        num_actual_inputs = len(tool._inputs)
        num_expected_inputs = 4
        self.assertEqual(num_actual_inputs, num_expected_inputs)
        
        tool = self.flat_wf.step_nodes['stp6'].tool
        
        # checking num inputs
        num_actual_inputs = len(tool._inputs)
        num_expected_inputs = 3
        self.assertEqual(num_actual_inputs, num_expected_inputs)
        
        # checking default values
        for tinput in tool._inputs:
            if tinput.id() == 'inStrOpt1':
                self.assertEqual(tinput.default, 'hello')
            elif tinput.id() == 'inStrOpt2':
                self.assertEqual(tinput.default, 'there')
            elif tinput.id() == 'inStrOpt3':
                self.assertEqual(tinput.default, 'friend')
            else:
                self.assertEqual(tinput.default, None)
        
    def test_remove_sources(self) -> None:
        step = self.flat_wf.step_nodes['stp4']
        actual_sources = set(list(step.sources.keys()))
        expected_sources = set(['inFileOpt1', 'inStrOpt1'])
        self.assertEqual(actual_sources, expected_sources)
        
        step = self.flat_wf.step_nodes['stp5']
        actual_sources = set(list(step.sources.keys()))
        expected_sources = set(['inFileOpt2', 'inStrOpt2'])
        self.assertEqual(actual_sources, expected_sources)
        
        step = self.flat_wf.step_nodes['stp6']
        actual_sources = set(list(step.sources.keys()))
        expected_sources = set()
        self.assertEqual(actual_sources, expected_sources)


    ### tool inputs: flat wf ###

    def test_mandatory_tinputs(self) -> None:
        tool = self.flat_wf.step_nodes['stp1'].tool
        actual_tinputs = set([x.id() for x in tool._inputs])
        expected_tinputs = set([
            'inFile1', 'inFile2', 'inStr1', 'inStr2'
        ])
        self.assertEqual(actual_tinputs, expected_tinputs)

    def test_connection_sources(self) -> None:
        tool = self.flat_wf.step_nodes['stp4'].tool
        actual_tinputs = set([x.id() for x in tool._inputs])
        expected_tinputs = set(['inFileOpt1', 'inStrOpt1'])
        for tinput_id in expected_tinputs:
            self.assertIn(tinput_id, actual_tinputs)
    
    def test_inputnode_sources(self) -> None:
        tool = self.flat_wf.step_nodes['stp5'].tool
        actual_tinputs = set([x.id() for x in tool._inputs])
        expected_tinputs = set(['inFileOpt2', 'inStrOpt2'])
        for tinput_id in expected_tinputs:
            self.assertIn(tinput_id, actual_tinputs)
    
    def test_static_sources(self) -> None:
        tool = self.flat_wf.step_nodes['stp6'].tool
        actual_tinputs = set([x.id() for x in tool._inputs])
        expected_tinputs = set(['inStrOpt3'])
        for tinput_id in expected_tinputs:
            self.assertIn(tinput_id, actual_tinputs)
    
    def test_optional_tinputs(self) -> None:
        tool = self.flat_wf.step_nodes['stp4'].tool
        actual_tinputs = set([x.id() for x in tool._inputs])
        expected_tinputs = set([
            'inFileOpt1', 
            'inFileOpt2', 
            'inStrOpt1', 
            'inStrOpt2', 
        ])
        self.assertEqual(actual_tinputs, expected_tinputs)
    
    def test_inputref_tinputs(self) -> None:
        tool = self.flat_wf.step_nodes['stp2'].tool
        actual_tinputs = set([x.id() for x in tool._inputs])
        expected_tinputs = set(['inFile1', 'OutputName'])
        self.assertEqual(actual_tinputs, expected_tinputs)
    
    def test_outputref_tinputs(self) -> None:
        tool = self.flat_wf.step_nodes['stp3'].tool
        actual_tinputs = set([x.id() for x in tool._inputs])
        expected_tinputs = set(['inFileOpt1', 'inStrOpt1'])
        self.assertEqual(actual_tinputs, expected_tinputs)
    
    ### tool inputs: nested wf ###
    
    def test_nested_tool(self) -> None:
        tool = self.nested_wf.step_nodes['stp1'].tool
        actual_tinputs = set([x.id() for x in tool._inputs])
        expected_tinputs = set([
            'inFileOpt1', 
            'inFileOpt2', 
            'inFileOpt3', 
            'inStrOpt1', 
            'inStrOpt2', 
            'inStrOpt3', 
        ])
        self.assertEqual(actual_tinputs, expected_tinputs)
    
    def test_nested_wf1(self) -> None:
        tool = self.nested_wf.step_nodes['stp2'].tool
        actual_tinputs = set([x for x in tool.input_nodes.keys()])
        expected_tinputs = set([
            'inFileOpt2', 
            'inFileOpt3', 
            'inStrOpt2', 
            'inStrOpt3', 
        ])
        self.assertEqual(actual_tinputs, expected_tinputs)

    def test_nested_wf2(self) -> None:
        tool = self.nested_wf.step_nodes['stp2'].tool.step_nodes['stp2'].tool
        actual_tinputs = set([x for x in tool.input_nodes.keys()])
        expected_tinputs = set([
            'inFileOpt3', 
            'inStrOpt3', 
        ])
        self.assertEqual(actual_tinputs, expected_tinputs)




# ---- PREPROCESSING: MODES ------------------------------

class TestPreprocessingModes(unittest.TestCase):
    
    def setUp(self) -> None:
        _reset_global_settings()

    def test_skeleton_nextflow(self) -> None:
        settings.translate.MODE = 'skeleton'
        filepath = f'{CWL_TESTDATA_PATH}/workflows/subworkflow_test/main.cwl'
        _, _, sub_tasks = _run(filepath, srcfmt='cwl', destfmt='nextflow')
        expected_inputs_count = {
            'modules/basic.nf': 3,
            'modules/mandatory_input_types.nf': 6,
            'modules/optional_input_types.nf': 5,
            'subworkflows/subworkflow.nf': 6,
        }
        expected_script_lengths = {
            'modules/basic.nf': 6,
            'modules/mandatory_input_types.nf': 8,
            'modules/optional_input_types.nf': 7,
        }
        for filepath, filecontents in sub_tasks:
            if _is_nf_process(filecontents):
                print(filecontents)
                actual_input_lines = _get_nf_process_input_lines(filecontents)
                actual_script_lines = _get_nf_process_script_lines(filecontents)
                self.assertEqual(len(actual_input_lines), expected_inputs_count[filepath])
                self.assertEqual(len(actual_script_lines), expected_script_lengths[filepath])

    def test_regular_nextflow(self) -> None:
        settings.translate.MODE = 'regular'
        filepath = f'{CWL_TESTDATA_PATH}/workflows/subworkflow_test/main.cwl'
        maintask, _, sub_tasks = _run(filepath, srcfmt='cwl', destfmt='nextflow')
        print(maintask)
        expected_inputs_count = {
            'modules/basic.nf': 3,
            'modules/mandatory_input_types.nf': 6,
            'modules/optional_input_types.nf': 5,
            'subworkflows/subworkflow.nf': 6,
        }
        expected_script_lengths = {
            'modules/basic.nf': 7,
            'modules/mandatory_input_types.nf': 8,
            'modules/optional_input_types.nf': 7,
        }
        for filepath, filecontents in sub_tasks:
            if _is_nf_process(filecontents):
                print(filecontents)
                actual_input_lines = _get_nf_process_input_lines(filecontents)
                actual_script_lines = _get_nf_process_script_lines(filecontents)
                self.assertEqual(len(actual_input_lines), expected_inputs_count[filepath])
                self.assertEqual(len(actual_script_lines), expected_script_lengths[filepath])

    def test_extended_nextflow(self) -> None:
        settings.translate.MODE = 'extended'
        filepath = f'{CWL_TESTDATA_PATH}/workflows/subworkflow_test/main.cwl'
        _, _, sub_tasks = _run(filepath, srcfmt='cwl', destfmt='nextflow')
        expected_inputs_count = {
            'modules/basic.nf': 7,
            'modules/mandatory_input_types.nf': 6,
            'modules/optional_input_types.nf': 6,
        }
        expected_script_lengths = {
            'modules/basic.nf': 10,
            'modules/mandatory_input_types.nf': 8,
            'modules/optional_input_types.nf': 8,
        }
        for filepath, filecontents in sub_tasks:
            if _is_nf_process(filecontents):
                print(filecontents)
                actual_input_lines = _get_nf_process_input_lines(filecontents)
                actual_script_lines = _get_nf_process_script_lines(filecontents)
                self.assertEqual(len(actual_input_lines), expected_inputs_count[filepath])
                self.assertEqual(len(actual_script_lines), expected_script_lengths[filepath])

    @unittest.skip('not implemented')
    def test_skeleton_cwl(self) -> None:
        settings.translate.MODE = 'skeleton'
        filepath = f'{CWL_TESTDATA_PATH}/workflows/subworkflow_test/main.cwl'
        maintask, _, sub_tasks = _run(filepath, srcfmt='cwl', destfmt='cwl')

        # main
        expected_num_clt_inputs = 11
        clt_inputs = _get_cwl_clt_inputs(maintask)
        self.assertEqual(len(clt_inputs), expected_num_clt_inputs)

        # subtasks
        expected_num_clt_inputs = {
            'tools/basic_v0_1_0.cwl': 4,
            'tools/mandatory_input_types_v0_1_0.cwl': 6,
            'tools/optional_input_types_v0_1_0.cwl': 5,
            'tools/subworkflow.cwl': 6,
        }
        for filepath, filecontents in sub_tasks:
            if _is_cwl_clt(filecontents):
                clt_inputs = _get_cwl_clt_inputs(filecontents)
                
                # checking expected number of clt inputs
                self.assertEqual(len(clt_inputs), expected_num_clt_inputs[filepath])

                # checking clt inputs have inputBindings
                for inp in clt_inputs:
                    self.assertNotIn('inputBinding', inp)
    
    @unittest.skip('not implemented')
    def test_skeleton_wdl(self) -> None:
        # TODO
        settings.translate.MODE = 'skeleton'
        filepath = f'{CWL_TESTDATA_PATH}/workflows/subworkflow_test/main.cwl'
        _, _, sub_tasks = _run(filepath, srcfmt='cwl', destfmt='wdl')
        for filepath, filecontents in sub_tasks:
            if _is_wdl_task(filecontents):
                command_lines = _get_wdl_task_command_lines(filecontents)
                self.assertEqual(len(command_lines), 2)
    
    @unittest.skip('not implemented')
    def test_regular_cwl1(self) -> None:
        settings.translate.MODE = 'regular'
        filepath = f'{CWL_TESTDATA_PATH}/workflows/subworkflow_test/main.cwl'
        maintask, _, sub_tasks = _run(filepath, srcfmt='cwl', destfmt='cwl')

        # main
        expected_num_clt_inputs = 12
        clt_inputs = _get_cwl_clt_inputs(maintask)
        self.assertEqual(len(clt_inputs), expected_num_clt_inputs)

        # subtasks
        expected_num_clt_inputs = {
            'tools/basic_v0_1_0.cwl': 6,
            'tools/mandatory_input_types_v0_1_0.cwl': 6,
            'tools/optional_input_types_v0_1_0.cwl': 5,
            'tools/subworkflow.cwl': 6,
        }
        for filepath, filecontents in sub_tasks:
            if _is_cwl_clt(filecontents):
                clt_inputs = _get_cwl_clt_inputs(filecontents)
                
                # checking expected number of clt inputs
                self.assertEqual(len(clt_inputs), expected_num_clt_inputs[filepath])

                # checking clt inputs have inputBindings
                for inp in clt_inputs:
                    self.assertIn('inputBinding', inp)
    
    @unittest.skip('not implemented')
    def test_regular_cwl2(self) -> None:
        settings.translate.MODE = 'regular'
        filepath = f'{CWL_TESTDATA_PATH}/workflows/m-unlock/workflows/ngtax.cwl'
        _, _, sub_tasks = _run(filepath, srcfmt='cwl', destfmt='cwl')
        expected_num_clt_inputs = {
            'tools/fastqc_v0_1_0.cwl': 2,
            'tools/files_to_folder_v0_1_0.cwl': 2,
            'tools/ngtax_v0_1_0.cwl': 9,
            'tools/ngtax_to_tsv_fasta_v0_1_0.cwl': 4,
        }
        expected_input_binding_absence = {
            'tools/fastqc_v0_1_0.cwl': [],
            'tools/files_to_folder_v0_1_0.cwl': ['files', 'folders', 'destination'],
            'tools/ngtax_v0_1_0.cwl': ['sample', 'fragment'],
            'tools/ngtax_to_tsv_fasta_v0_1_0.cwl': ['metadata', 'input', 'identifier', 'fragment'],
        }
        expected_num_clt_args = {
            'tools/fastqc_v0_1_0.cwl': 2,
            'tools/files_to_folder_v0_1_0.cwl': 1,
            'tools/ngtax_v0_1_0.cwl': 4,
            'tools/ngtax_to_tsv_fasta_v0_1_0.cwl': 8,
        }
        for filepath, filecontents in sub_tasks:
            if _is_cwl_clt(filecontents):
                clt_inputs = _get_cwl_clt_inputs(filecontents)
                clt_args = _get_cwl_clt_args(filecontents)
                
                # checking expected number of clt inputs
                self.assertEqual(len(clt_inputs), expected_num_clt_inputs[filepath])
                # checking expected number of clt args
                self.assertEqual(len(clt_args), expected_num_clt_args[filepath])

                # checking clt inputs have or do not have inputBindings
                for inp in clt_inputs:
                    if inp['id'] not in expected_input_binding_absence[filepath]:
                        self.assertIn('inputBinding', inp)
                    else:
                        self.assertNotIn('inputBinding', inp)
    
    @unittest.skip('implement me')
    def test_regular_wdl(self) -> None:
        settings.translate.MODE = 'regular'
        filepath = f'{CWL_TESTDATA_PATH}/workflows/subworkflow_test/main.cwl'
        _, _, sub_tasks = _run(filepath, srcfmt='cwl', destfmt='wdl')
        expected_num_clt_inputs = {
            'align_and_tag_v0_1_0': 3,
            'index_bam_v0_1_0': 1,
            'mark_duplicates_and_sort_v0_1_0': 3,
            'merge_bams_samtools_v0_1_0': 2,
            'name_sort_v0_1_0': 1,
        }
        expected_input_binding_absence = {
            'align_and_tag_v0_1_0': [],
            'index_bam_v0_1_0': ['bam'],
            'mark_duplicates_and_sort_v0_1_0': [],
            'merge_bams_samtools_v0_1_0': ['name'],
            'name_sort_v0_1_0': [],
        }
        for filepath, filecontents in sub_tasks:
            if _is_wdl_task(filecontents):
                command_lines = _get_wdl_task_command_lines(filecontents)
                raise NotImplementedError

    @unittest.skip('not implemented')
    def test_extended_cwl(self) -> None:
        settings.translate.MODE = 'extended'
        filepath = f'{CWL_TESTDATA_PATH}/workflows/subworkflow_test/main.cwl'
        _, _, sub_tasks = _run(filepath, srcfmt='cwl', destfmt='cwl')
        expected_num_clt_inputs = {
            'tools/basic_v0_1_0.cwl': 7,
            'tools/mandatory_input_types_v0_1_0.cwl': 6,
            'tools/optional_input_types_v0_1_0.cwl': 6,
        }
        for filepath, filecontents in sub_tasks:
            if _is_cwl_clt(filecontents):
                clt_inputs = _get_cwl_clt_inputs(filecontents)
                
                # checking expected number of clt inputs
                self.assertEqual(len(clt_inputs), expected_num_clt_inputs[filepath])

                # checking clt inputs have inputBindings
                for inp in clt_inputs:
                    self.assertIn('inputBinding', inp)

    @unittest.skip('implement me')
    def test_extended_wdl(self) -> None:
        settings.translate.MODE = 'extended'
        filepath = f'{CWL_TESTDATA_PATH}/workflows/subworkflow_test/main.cwl'
        _, _, sub_tasks = _run(filepath, srcfmt='cwl', destfmt='wdl')
        expected_num_clt_inputs = {
            'align_and_tag_v0_1_0': 3,
            'index_bam_v0_1_0': 1,
            'mark_duplicates_and_sort_v0_1_0': 3,
            'merge_bams_samtools_v0_1_0': 2,
            'name_sort_v0_1_0': 1,
        }
        expected_input_binding_absence = {
            'align_and_tag_v0_1_0': [],
            'index_bam_v0_1_0': ['bam'],
            'mark_duplicates_and_sort_v0_1_0': [],
            'merge_bams_samtools_v0_1_0': ['name'],
            'name_sort_v0_1_0': [],
        }
        for filepath, filecontents in sub_tasks:
            if _is_wdl_task(filecontents):
                raise NotImplementedError



# ---- PREPROCESSING: TO BUILDERS ------------------------------



class TestPreprocessingToBuilders(unittest.TestCase):
    
    def setUp(self) -> None:
        _reset_global_settings()

    def test_commandtool_to_builders(self) -> None:
        entity = GridssTestTool()
        entity = to_builders(entity)
        self.assertIsInstance(entity, CommandToolBuilder)
    
    def test_codetool_to_builders(self) -> None:
        entity = FileOutputPythonTestTool()
        entity = to_builders(entity)
        self.assertIsInstance(entity, CodeTool)

    def test_workflow_to_builders(self) -> None:
        entity = AssemblyTestWF()
        entity = to_builders(entity)
        self.assertIsInstance(entity, WorkflowBuilder)
        assert(isinstance(entity, WorkflowBuilder))
        for step in entity.step_nodes.values():
            self.assertIsInstance(step.tool, CommandToolBuilder)



# ---- FROM CWL ---------------------------

class TestCwlToWdl(unittest.TestCase):
    
    def setUp(self) -> None:
        self.src = 'cwl'
        self.dest = 'wdl'
        _reset_global_settings()

    def test_tool_fastqc2(self):
        filepath = f'{CWL_TESTDATA_PATH}/tools/fastqc2.cwl'
        toolstr = _run(filepath, self.src, self.dest)
        print(toolstr)
    
    def test_wf_super_enhancer(self):
        filepath = f'{CWL_TESTDATA_PATH}/workflows/super_enhancer_wf.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)

    @unittest.skip('implement secondary type mismatch cleanup')
    def test_wf_kids_manta(self):
        filepath = f'{CWL_TESTDATA_PATH}/workflows/kf-somatic-workflow/workflow/kfdrc_production_manta_wf.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    @unittest.skip('implement scatter on multiple inputs')
    def test_wf_ebi_metagenomics_raw_reads(self):
        filepath = f'{CWL_TESTDATA_PATH}/workflows/ebi-metagenomics/workflows/raw-reads-wf--v.5-cond.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    @unittest.skip('implement scatter on multiple inputs')
    def test_wf_ebi_metagenomics_amplicon(self):
        filepath = f'{CWL_TESTDATA_PATH}/workflows/ebi-metagenomics/workflows/amplicon-wf--v.5-cond.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    @unittest.skip('implement scatter on multiple inputs')
    def test_wf_ebi_metagenomics_assembly(self):
        filepath = f'{CWL_TESTDATA_PATH}/workflows/ebi-metagenomics/workflows/assembly-wf--v.5-cond.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)

    def test_wf_munlock_demultiplexing(self):
        filepath = f'{CWL_TESTDATA_PATH}/workflows/m-unlock/workflows/demultiplexing.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    def test_wf_munlock_mock_ngtax(self):
        filepath = f'{CWL_TESTDATA_PATH}/workflows/m-unlock/workflows/mock_ngtax.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    def test_wf_munlock_pilon_mapping(self):
        filepath = f'{CWL_TESTDATA_PATH}/workflows/m-unlock/workflows/pilon_mapping.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    def test_wf_munlock_sapp_microbes(self):
        filepath = f'{CWL_TESTDATA_PATH}/workflows/m-unlock/workflows/sapp_microbes.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)

    def test_wf_munlock_toHDT_compression(self):
        filepath = f'{CWL_TESTDATA_PATH}/workflows/m-unlock/workflows/toHDT_compression.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    def test_wf_munlock_ngtax(self):
        filepath = f'{CWL_TESTDATA_PATH}/workflows/m-unlock/workflows/ngtax.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
        
    def test_wf_munlock_metagenomics_GEM(self):
        filepath = f'{CWL_TESTDATA_PATH}/workflows/m-unlock/workflows/metagenomics_GEM.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)

    def test_wf_munlock_ngtax_picrust2(self):
        filepath = f'{CWL_TESTDATA_PATH}/workflows/m-unlock/workflows/ngtax_picrust2.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)

    # hard
    def test_wf_cromast(self):
        filepath = f'{CWL_TESTDATA_PATH}/workflows/CroMaSt/CroMaSt.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)


class TestCwlToNextflow(unittest.TestCase):
    
    def setUp(self) -> None:
        self.src = 'cwl'
        self.dest = 'nextflow'
        _reset_global_settings()

    # Tools
    def test_tool_samtools_flagstat(self):
        filepath = f'{CWL_TESTDATA_PATH}/tools/samtools_flagstat.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)

    def test_tool_gatk_haplotype_caller(self):
        filepath = f'{CWL_TESTDATA_PATH}/tools/gatk_haplotype_caller.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)

    def test_tool_fastqc(self):
        filepath = f'{CWL_TESTDATA_PATH}/tools/fastqc2.cwl'
        toolstr = _run(filepath, self.src, self.dest)
        settings.translate.nextflow.ENTITY = 'workflow'
        print(toolstr)
    
    # Workflows
    def test_wf_align_sort_markdup(self):
        filepath = f'{CWL_TESTDATA_PATH}/workflows/align_sort_markdup/align_sort_markdup.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)

    def test_wf_super_enhancer(self):
        filepath = f'{CWL_TESTDATA_PATH}/workflows/super_enhancer_wf.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)

    def test_wf_kids_manta(self):
        filepath = f'{CWL_TESTDATA_PATH}/workflows/kf-somatic-workflow/workflow/kfdrc_production_manta_wf.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
     
    def test_wf_ebi_metagenomics_raw_reads(self):
        filepath = f'{CWL_TESTDATA_PATH}/workflows/ebi-metagenomics/workflows/raw-reads-wf--v.5-cond.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    def test_wf_ebi_metagenomics_amplicon(self):
        filepath = f'{CWL_TESTDATA_PATH}/workflows/ebi-metagenomics/workflows/amplicon-wf--v.5-cond.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    def test_wf_ebi_metagenomics_assembly(self):
        filepath = f'{CWL_TESTDATA_PATH}/workflows/ebi-metagenomics/workflows/assembly-wf--v.5-cond.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)

    def test_wf_munlock_demultiplexing(self):
        filepath = f'{CWL_TESTDATA_PATH}/workflows/m-unlock/workflows/demultiplexing.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    def test_wf_munlock_mock_ngtax(self):
        filepath = f'{CWL_TESTDATA_PATH}/workflows/m-unlock/workflows/mock_ngtax.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    def test_wf_munlock_pilon_mapping(self):
        filepath = f'{CWL_TESTDATA_PATH}/workflows/m-unlock/workflows/pilon_mapping.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    def test_wf_munlock_sapp_microbes(self):
        filepath = f'{CWL_TESTDATA_PATH}/workflows/m-unlock/workflows/sapp_microbes.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)

    def test_wf_munlock_toHDT_compression(self):
        filepath = f'{CWL_TESTDATA_PATH}/workflows/m-unlock/workflows/toHDT_compression.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    def test_wf_munlock_ngtax(self):
        filepath = f'{CWL_TESTDATA_PATH}/workflows/m-unlock/workflows/ngtax.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
        
    def test_wf_munlock_metagenomics_GEM(self):
        filepath = f'{CWL_TESTDATA_PATH}/workflows/m-unlock/workflows/metagenomics_GEM.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)

    def test_wf_munlock_ngtax_picrust2(self):
        filepath = f'{CWL_TESTDATA_PATH}/workflows/m-unlock/workflows/ngtax_picrust2.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)

    # hard
    @unittest.skip('need to check prune to make sure cwl entities have unique names')
    def test_wf_cromast(self):
        filepath = f'{CWL_TESTDATA_PATH}/workflows/CroMaSt/CroMaSt.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)



# ---- FROM WDL ---------------------------

class TestWdlToCwl(unittest.TestCase):
    
    def setUp(self) -> None:
        self.src = 'wdl'
        self.dest = 'cwl'
        _reset_global_settings()

    @unittest.skip('need injest fixes')
    def test_tool_bwa(self):
        filepath = f'{WDL_TESTDATA_PATH}/bwa.wdl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    @unittest.skip('need injest fixes')
    def test_wf_somatic(self):
        filepath = f'{WDL_TESTDATA_PATH}/somatic_wf.wdl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    @unittest.skip('need injest fixes')
    def test_wf_multisample_jointgt_gatk4_wf(self):
        filepath = f'{WDL_TESTDATA_PATH}/Multisample_jointgt_GATK4.wdl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)

    @unittest.skip('need injest fixes')
    def test_wf_reads2map_preprocessing(self):
        filepath = f'{WDL_TESTDATA_PATH}/Reads2Map/pipelines/PreprocessingReads/PreprocessingReads.wdl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)

    @unittest.skip('need injest fixes')
    def test_wf_reads2map_reads2map(self):
        filepath = f'{WDL_TESTDATA_PATH}/Reads2Map/pipelines/EmpiricalReads2Map/EmpiricalReads2Map.wdl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    @unittest.skip('need injest fixes')
    def test_wf_reads2map_snp_calling(self):
        filepath = f'{WDL_TESTDATA_PATH}/Reads2Map/pipelines/EmpiricalSNPCalling/EmpiricalSNPCalling.wdl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)

    


class TestWdlToNextflow(unittest.TestCase):
    
    def setUp(self) -> None:
        self.src = 'wdl'
        self.dest = 'nextflow'
        _reset_global_settings()

    @unittest.skip('need injest fixes')
    def test_wf_multisample_jointgt_gatk4(self):
        filepath = f'{WDL_TESTDATA_PATH}/Multisample_jointgt_GATK4.wdl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)

    @unittest.skip('need injest fixes')
    def test_wf_reads2map_preprocessing(self):
        filepath = f'{WDL_TESTDATA_PATH}/Reads2Map/pipelines/PreprocessingReads/PreprocessingReads.wdl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)

    @unittest.skip('need injest fixes')
    def test_wf_reads2map_reads2map(self):
        filepath = f'{WDL_TESTDATA_PATH}/Reads2Map/pipelines/EmpiricalReads2Map/EmpiricalReads2Map.wdl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    @unittest.skip('need injest fixes')
    def test_wf_reads2map_snp_calling(self):
        filepath = f'{WDL_TESTDATA_PATH}/Reads2Map/pipelines/EmpiricalSNPCalling/EmpiricalSNPCalling.wdl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)



# ---- FROM GALAXY ------------------------

class TestGalaxyToWdl(unittest.TestCase):
    
    def setUp(self) -> None:
        self.src = 'galaxy'
        self.dest = 'wdl'
        _reset_global_settings()
        settings.translate.MODE = 'regular'
    
    def test_wf_abricate(self):
        filepath = f'{GALAXY_TESTDATA_PATH}/wf_abricate.ga'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)

    def test_wf_unicycler_assembly(self):
        filepath = f'{GALAXY_TESTDATA_PATH}/unicycler_assembly.ga'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)

    def test_wf_rna_seq_counts_to_genes(self):
        filepath = f'{GALAXY_TESTDATA_PATH}/rna_seq_counts_to_genes.ga'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    def test_wf_rna_seq_genes_to_pathways(self):
        filepath = f'{GALAXY_TESTDATA_PATH}/rna_seq_genes_to_pathways.ga'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    @unittest.skip('implement scatter on multiple inputs')
    def test_wf_rna_seq_reads_to_counts(self):
        filepath = f'{GALAXY_TESTDATA_PATH}/rna_seq_reads_to_counts.ga'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)



class TestGalaxyToCwl(unittest.TestCase):
    
    def setUp(self) -> None:
        self.src = 'galaxy'
        self.dest = 'cwl'
        _reset_global_settings()
        settings.translate.MODE = 'regular'

    def test_wf_abricate(self):
        filepath = f'{GALAXY_TESTDATA_PATH}/wf_abricate.ga'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)

    def test_wf_unicycler_assembly(self):
        filepath = f'{GALAXY_TESTDATA_PATH}/unicycler_assembly.ga'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)

    def test_wf_rna_seq_counts_to_genes(self):
        filepath = f'{GALAXY_TESTDATA_PATH}/rna_seq_counts_to_genes.ga'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    def test_wf_rna_seq_genes_to_pathways(self):
        filepath = f'{GALAXY_TESTDATA_PATH}/rna_seq_genes_to_pathways.ga'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    def test_wf_rna_seq_reads_to_counts(self):
        filepath = f'{GALAXY_TESTDATA_PATH}/rna_seq_reads_to_counts.ga'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)


