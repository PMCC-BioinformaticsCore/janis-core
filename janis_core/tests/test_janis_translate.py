




import unittest
from typing import Optional

from janis_core.ingestion import ingest
from janis_core.translations import translate

from janis_core import settings
# PATHS MUST BE ABSOLUTE




# ------- HELPER FUNCS ------- #

def _run(filepath: str, srcfmt: str, destfmt: str) -> Optional[str]:
    internal = ingest(filepath, srcfmt)
    return translate(internal, destfmt, allow_empty_container=True, export_path='./translated')

def _reset_global_settings() -> None:
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
        filepath = '/home/grace/work/pp/translation/examples/analysis-workflows/definitions/tools/samtools_flagstat.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    def test_tool_gatk_haplotype_caller(self):
        filepath = '/home/grace/work/pp/translation/examples/analysis-workflows/definitions/tools/gatk_haplotype_caller.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    def test_wf_align_sort_markdup(self):
        filepath = '/home/grace/work/pp/translation/examples/analysis-workflows/definitions/subworkflows/align_sort_markdup.cwl'
        mainstr = _run(filepath, self.src, self.dest)
        print(mainstr)
    
    def test_wf_alignment_exome(self):
        filepath = '/home/grace/work/pp/translation/examples/analysis-workflows/definitions/pipelines/alignment_exome.cwl'
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


