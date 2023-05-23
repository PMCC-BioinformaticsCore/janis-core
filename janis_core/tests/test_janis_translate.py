




import unittest
from typing import Optional, Any, Tuple

from janis_core.ingestion import ingest
from janis_core.translations import translate

from janis_core import settings
# PATHS MUST BE ABSOLUTE

import regex as re
import yaml


# ------- HELPER FUNCS ------- #

def _run(filepath: str, srcfmt: str, destfmt: str) -> Any:
    internal = ingest(filepath, srcfmt)
    return translate(internal, destfmt, allow_empty_container=True, export_path='./translated')

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
    settings.translate.MODE = 'minimal'
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


# ------- WORSHOP TESTS ------- #

class TestWorkshopCwlToNextflow(unittest.TestCase):
    
    def setUp(self) -> None:
        self.src = 'cwl'
        self.dest = 'nextflow'
        _reset_global_settings()

    def test_tool_samtools_flagstat(self):
        filepath = './janis_core/tests/data/cwl/workflows/analysis-workflows/tools/samtools_flagstat.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    def test_tool_gatk_haplotype_caller(self):
        filepath = './janis_core/tests/data/cwl/workflows/analysis-workflows/tools/gatk_haplotype_caller.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    def test_wf_align_sort_markdup(self):
        settings.translate.MODE = 'full'
        filepath = './janis_core/tests/data/cwl/workflows/analysis-workflows/subworkflows/align_sort_markdup.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    def test_wf_alignment_exome(self):
        filepath = './janis_core/tests/data/cwl/workflows/analysis-workflows/pipelines/alignment_exome.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)


class TestWorkshopGalaxyToNextflow(unittest.TestCase):
    
    def setUp(self) -> None:
        self.src = 'galaxy'
        self.dest = 'nextflow'
        _reset_global_settings()

    def test_abricate_wf(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/galaxy/wf_abricate.ga'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)

    def test_unicycler_assembly(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/galaxy/unicycler_assembly.ga'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)

    def test_rna_seq_counts_to_genes(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/galaxy/rna_seq_counts_to_genes.ga'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    def test_rna_seq_genes_to_pathways(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/galaxy/rna_seq_genes_to_pathways.ga'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    def test_rna_seq_reads_to_counts(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/galaxy/rna_seq_reads_to_counts.ga'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)


# ---- MODES ------------------------------

class TestModes(unittest.TestCase):
    
    def setUp(self) -> None:
        _reset_global_settings()
    
    def test_skeleton_cwl(self) -> None:
        settings.translate.MODE = 'skeleton'
        filepath = './janis_core/tests/data/cwl/workflows/subworkflow_test/main.cwl'
        _, _, sub_tasks = _run(filepath, srcfmt='cwl', destfmt='cwl')
        expected_num_clt_inputs = {
            'tools/basic_v0_1_0.cwl': 4,
            'tools/mandatory_input_types_v0_1_0.cwl': 6,
            'tools/optional_input_types_v0_1_0.cwl': 5,
        }
        for filepath, filecontents in sub_tasks:
            if _is_cwl_clt(filecontents):
                clt_inputs = _get_cwl_clt_inputs(filecontents)
                
                # checking expected number of clt inputs
                self.assertEqual(len(clt_inputs), expected_num_clt_inputs[filepath])

                # checking clt inputs have inputBindings
                for inp in clt_inputs:
                    self.assertNotIn('inputBinding', inp)
    
    def test_skeleton_wdl(self) -> None:
        settings.translate.MODE = 'skeleton'
        filepath = './janis_core/tests/data/cwl/workflows/subworkflow_test/main.cwl'
        _, _, sub_tasks = _run(filepath, srcfmt='cwl', destfmt='wdl')
        for filepath, filecontents in sub_tasks:
            if _is_wdl_task(filecontents):
                command_lines = _get_wdl_task_command_lines(filecontents)
                self.assertEqual(len(command_lines), 2)
    
    def test_skeleton_nextflow(self) -> None:
        settings.translate.MODE = 'skeleton'
        filepath = './janis_core/tests/data/cwl/workflows/subworkflow_test/main.cwl'
        _, _, sub_tasks = _run(filepath, srcfmt='cwl', destfmt='nextflow')
        expected_inputs_count = {
            'modules/basic.nf': 4,
            'modules/mandatory_input_types.nf': 6,
            'modules/optional_input_types.nf': 5,
        }
        expected_script_lengths = {
            'modules/basic.nf': 1,
            'modules/mandatory_input_types.nf': 1,
            'modules/optional_input_types.nf': 1,
        }
        for filepath, filecontents in sub_tasks:
            if _is_nf_process(filecontents):
                actual_input_lines = _get_nf_process_input_lines(filecontents)
                actual_script_lines = _get_nf_process_script_lines(filecontents)
                self.assertEqual(len(actual_input_lines), expected_inputs_count[filepath])
                self.assertEqual(len(actual_script_lines), expected_script_lengths[filepath])
    
    def test_minimal_cwl(self) -> None:
        settings.translate.MODE = 'minimal'
        filepath = './janis_core/tests/data/cwl/workflows/subworkflow_test/main.cwl'
        _, _, sub_tasks = _run(filepath, srcfmt='cwl', destfmt='cwl')
        expected_num_clt_inputs = {
            'tools/basic_v0_1_0.cwl': 4,
            'tools/mandatory_input_types_v0_1_0.cwl': 6,
            'tools/optional_input_types_v0_1_0.cwl': 5,
        }
        for filepath, filecontents in sub_tasks:
            if _is_cwl_clt(filecontents):
                clt_inputs = _get_cwl_clt_inputs(filecontents)
                
                # checking expected number of clt inputs
                self.assertEqual(len(clt_inputs), expected_num_clt_inputs[filepath])

                # checking clt inputs have inputBindings
                for inp in clt_inputs:
                    self.assertIn('inputBinding', inp)
    
    def test_minimal_wdl(self) -> None:
        settings.translate.MODE = 'minimal'
        filepath = './janis_core/tests/data/cwl/workflows/subworkflow_test/main.cwl'
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
            if _is_cwl_clt(filecontents):
                command_lines = _get_wdl_task_command_lines(filecontents)
                
                # checking expected number of clt inputs
                toolname = filepath.split('/')[-1].split('.')[0]
                self.assertEqual(len(clt_inputs), expected_num_clt_inputs[toolname])

                # checking clt inputs have inputBindings
                for inp in clt_inputs:
                    if inp['id'] in expected_input_binding_absence[toolname]:
                        continue
                    else:
                        self.assertIn('inputBinding', inp)
    
    def test_minimal_nextflow(self) -> None:
        settings.translate.MODE = 'minimal'
        filepath = './janis_core/tests/data/cwl/workflows/subworkflow_test/main.cwl'
        _, _, sub_tasks = _run(filepath, srcfmt='cwl', destfmt='nextflow')
        expected_inputs_count = {
            'modules/basic.nf': 4,
            'modules/mandatory_input_types.nf': 6,
            'modules/optional_input_types.nf': 5,
        }
        expected_script_lengths = {
            'modules/basic.nf': 6,
            'modules/mandatory_input_types.nf': 7,
            'modules/optional_input_types.nf': 6,
        }
        for filepath, filecontents in sub_tasks:
            if _is_nf_process(filecontents):
                actual_input_lines = _get_nf_process_input_lines(filecontents)
                actual_script_lines = _get_nf_process_script_lines(filecontents)
                self.assertEqual(len(actual_input_lines), expected_inputs_count[filepath])
                self.assertEqual(len(actual_script_lines), expected_script_lengths[filepath])

    def test_full_cwl(self) -> None:
        settings.translate.MODE = 'full'
        filepath = './janis_core/tests/data/cwl/workflows/subworkflow_test/main.cwl'
        _, _, sub_tasks = _run(filepath, srcfmt='cwl', destfmt='cwl')
        expected_num_clt_inputs = {
            'tools/basic_v0_1_0.cwl': 6,
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

    def test_full_nextflow(self) -> None:
        settings.translate.MODE = 'full'
        filepath = './janis_core/tests/data/cwl/workflows/subworkflow_test/main.cwl'
        _, _, sub_tasks = _run(filepath, srcfmt='cwl', destfmt='nextflow')
        expected_inputs_count = {
            'modules/basic.nf': 6,
            'modules/mandatory_input_types.nf': 6,
            'modules/optional_input_types.nf': 6,
        }
        expected_script_lengths = {
            'modules/basic.nf': 8,
            'modules/mandatory_input_types.nf': 7,
            'modules/optional_input_types.nf': 7,
        }
        for filepath, filecontents in sub_tasks:
            if _is_nf_process(filecontents):
                actual_input_lines = _get_nf_process_input_lines(filecontents)
                actual_script_lines = _get_nf_process_script_lines(filecontents)
                self.assertEqual(len(actual_input_lines), expected_inputs_count[filepath])
                self.assertEqual(len(actual_script_lines), expected_script_lengths[filepath])



# ---- FROM CWL ---------------------------

class TestCwlToWdl(unittest.TestCase):
    
    def setUp(self) -> None:
        self.src = 'cwl'
        self.dest = 'wdl'
        _reset_global_settings()
    
    def test_super_enhancer(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/super_enhancer_wf.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)

    def test_fastqc2_tool(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/fastqc2.cwl'
        toolstr = _run(filepath, self.src, self.dest)
        print(toolstr)

    def test_kids_manta(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/kf-somatic-workflow/workflow/kfdrc_production_manta_wf.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    def test_ebi_metagenomics_raw_reads(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/ebi-metagenomics/workflows/raw-reads-wf--v.5-cond.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    def test_ebi_metagenomics_amplicon(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/ebi-metagenomics/workflows/amplicon-wf--v.5-cond.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    def test_ebi_metagenomics_assembly(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/ebi-metagenomics/workflows/assembly-wf--v.5-cond.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)

    def test_munlock_demultiplexing(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/m-unlock/workflows/demultiplexing.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    def test_munlock_mock_ngtax(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/m-unlock/workflows/mock_ngtax.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    def test_munlock_pilon_mapping(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/m-unlock/workflows/pilon_mapping.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    def test_munlock_sapp_microbes(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/m-unlock/workflows/sapp_microbes.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)

    def test_munlock_toHDT_compression(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/m-unlock/workflows/toHDT_compression.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    def test_munlock_ngtax(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/m-unlock/workflows/ngtax.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
        
    def test_munlock_metagenomics_GEM(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/m-unlock/workflows/metagenomics_GEM.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)

    def test_munlock_ngtax_picrust2(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/m-unlock/workflows/ngtax_picrust2.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)

    # hard
    def test_cromast(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/CroMaSt/CroMaSt.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)


class TestCwlToNextflow(unittest.TestCase):
    
    def setUp(self) -> None:
        self.src = 'cwl'
        self.dest = 'nextflow'
        _reset_global_settings()

    def test_super_enhancer(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/super_enhancer_wf.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    def test_fastqc2_tool(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/fastqc2.cwl'
        toolstr = _run(filepath, self.src, self.dest)
        from janis_core import settings
        settings.translate.nextflow.MODE = 'workflow'
        print(toolstr)

    def test_kids_manta(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/kf-somatic-workflow/workflow/kfdrc_production_manta_wf.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
     
    def test_ebi_metagenomics_raw_reads(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/ebi-metagenomics/workflows/raw-reads-wf--v.5-cond.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    def test_ebi_metagenomics_amplicon(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/ebi-metagenomics/workflows/amplicon-wf--v.5-cond.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    def test_ebi_metagenomics_assembly(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/ebi-metagenomics/workflows/assembly-wf--v.5-cond.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)

    def test_munlock_demultiplexing(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/m-unlock/workflows/demultiplexing.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    def test_munlock_mock_ngtax(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/m-unlock/workflows/mock_ngtax.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    def test_munlock_pilon_mapping(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/m-unlock/workflows/pilon_mapping.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    def test_munlock_sapp_microbes(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/m-unlock/workflows/sapp_microbes.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)

    def test_munlock_toHDT_compression(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/m-unlock/workflows/toHDT_compression.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    def test_munlock_ngtax(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/m-unlock/workflows/ngtax.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
        
    def test_munlock_metagenomics_GEM(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/m-unlock/workflows/metagenomics_GEM.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)

    def test_munlock_ngtax_picrust2(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/m-unlock/workflows/ngtax_picrust2.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)

    # hard
    def test_cromast(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/CroMaSt/CroMaSt.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)



# ---- FROM WDL ---------------------------

class TestWdlToCwl(unittest.TestCase):
    
    def setUp(self) -> None:
        self.src = 'wdl'
        self.dest = 'cwl'
        _reset_global_settings()

    def test_multisample_jointgt_gatk4(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/wdl/Multisample_jointgt_GATK4.wdl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)

    def test_reads2map_preprocessing(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/wdl/Reads2Map/pipelines/PreprocessingReads/PreprocessingReads.wdl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)

    def test_reads2map_reads2map(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/wdl/Reads2Map/pipelines/EmpiricalReads2Map/EmpiricalReads2Map.wdl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    def test_reads2map_snp_calling(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/wdl/Reads2Map/pipelines/EmpiricalSNPCalling/EmpiricalSNPCalling.wdl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)

    


class TestWdlToNextflow(unittest.TestCase):
    
    def setUp(self) -> None:
        self.src = 'wdl'
        self.dest = 'nextflow'
        _reset_global_settings()

    def test_multisample_jointgt_gatk4(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/wdl/Multisample_jointgt_GATK4.wdl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)

    def test_reads2map_preprocessing(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/wdl/Reads2Map/pipelines/PreprocessingReads/PreprocessingReads.wdl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)

    def test_reads2map_reads2map(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/wdl/Reads2Map/pipelines/EmpiricalReads2Map/EmpiricalReads2Map.wdl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    def test_reads2map_snp_calling(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/wdl/Reads2Map/pipelines/EmpiricalSNPCalling/EmpiricalSNPCalling.wdl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)





# ---- FROM GALAXY ------------------------

class TestGalaxyToWdl(unittest.TestCase):
    
    def setUp(self) -> None:
        self.src = 'galaxy'
        self.dest = 'wdl'
        _reset_global_settings()
    
    def test_abricate_wf(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/galaxy/wf_abricate.ga'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)

    def test_unicycler_assembly(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/galaxy/unicycler_assembly.ga'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)

    def test_rna_seq_counts_to_genes(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/galaxy/rna_seq_counts_to_genes.ga'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    def test_rna_seq_genes_to_pathways(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/galaxy/rna_seq_genes_to_pathways.ga'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    def test_rna_seq_reads_to_counts(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/galaxy/rna_seq_reads_to_counts.ga'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)



class TestGalaxyToCwl(unittest.TestCase):
    
    def setUp(self) -> None:
        self.src = 'galaxy'
        self.dest = 'cwl'
        _reset_global_settings()

    def test_abricate_wf(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/galaxy/wf_abricate.ga'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)

    def test_unicycler_assembly(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/galaxy/unicycler_assembly.ga'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)

    def test_rna_seq_counts_to_genes(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/galaxy/rna_seq_counts_to_genes.ga'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    def test_rna_seq_genes_to_pathways(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/galaxy/rna_seq_genes_to_pathways.ga'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    def test_rna_seq_reads_to_counts(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/galaxy/rna_seq_reads_to_counts.ga'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)


