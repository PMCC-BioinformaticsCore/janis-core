
import unittest

from janis_core.ingestion import ingest
from janis_core.translations import NextflowTranslator


class TestCwlToNextflow(unittest.TestCase):
    
    def setUp(self) -> None:
        pass

    def test_munlock_demultiplexing(self):
        filepath = 'janis_core/tests/data/cwl/m-unlock/workflows/workflow_demultiplexing.cwl'
        wf = ingest(filepath, 'cwl')
        mainstr, substr_dict = NextflowTranslator.translate_workflow(wf)
        print('\n\nmain workflow -----------------\n')
        print(mainstr)
        for path, contents in substr_dict.items():
            print(f'\n\n{path} -----------------\n')
            print(contents)
        print()
        
    def test_munlock_metagenomics_GEM(self):
        filepath = 'janis_core/tests/data/cwl/m-unlock/workflows/workflow_metagenomics_GEM.cwl'
        wf = ingest(filepath, 'cwl')
        mainstr, substr_dict = NextflowTranslator.translate_workflow(wf)
        
    def test_munlock_mock_ngtax(self):
        filepath = 'janis_core/tests/data/cwl/m-unlock/workflows/workflow_mock_ngtax.cwl'
        wf = ingest(filepath, 'cwl')
        mainstr, substr_dict = NextflowTranslator.translate_workflow(wf)
        
    def test_munlock_ngtax_picrust2(self):
        filepath = 'janis_core/tests/data/cwl/m-unlock/workflows/workflow_ngtax_picrust2.cwl'
        wf = ingest(filepath, 'cwl')
        mainstr, substr_dict = NextflowTranslator.translate_workflow(wf)
        
    def test_munlock_ngtax(self):
        filepath = 'janis_core/tests/data/cwl/m-unlock/workflows/workflow_ngtax.cwl'
        wf = ingest(filepath, 'cwl')
        mainstr, substr_dict = NextflowTranslator.translate_workflow(wf)
        
    def test_munlock_pilon_mapping(self):
        filepath = 'janis_core/tests/data/cwl/m-unlock/workflows/workflow_pilon_mapping.cwl'
        wf = ingest(filepath, 'cwl')
        mainstr, substr_dict = NextflowTranslator.translate_workflow(wf)
        
    def test_munlock_sapp_microbes(self):
        filepath = 'janis_core/tests/data/cwl/m-unlock/workflows/workflow_sapp_microbes.cwl'
        wf = ingest(filepath, 'cwl')
        mainstr, substr_dict = NextflowTranslator.translate_workflow(wf)
        
    def test_munlock_toHDT_compression(self):
        filepath = 'janis_core/tests/data/cwl/m-unlock/workflows/workflow_toHDT_compression.cwl'
        wf = ingest(filepath, 'cwl')
        mainstr, substr_dict = NextflowTranslator.translate_workflow(wf)
    
    def test_cromast(self):
        filepath = 'janis_core/tests/data/cwl/CroMaSt/CroMaSt.cwl'
        wf = ingest(filepath, 'cwl')
        mainstr, substr_dict = NextflowTranslator.translate_workflow(wf)
    
    # MetaGOflow extremely broken
    def test_metagoflow(self):
        filepath = 'janis_core/tests/data/cwl/MetaGOflow/workflows/gos_wf.cwl'
        wf = ingest(filepath, 'cwl')
        mainstr, substr_dict = NextflowTranslator.translate_workflow(wf)
    
    def test_metagoflow_conditionals_qc(self):
        filepath = 'janis_core/tests/data/cwl/MetaGOflow/workflows/conditionals/qc.cwl'
        wf = ingest(filepath, 'cwl')
        mainstr, substr_dict = NextflowTranslator.translate_workflow(wf)
    
    def test_metagoflow_conditionals_rawreads(self):
        filepath = 'janis_core/tests/data/cwl/MetaGOflow/workflows/conditionals/raw-reads-2.cwl'
        wf = ingest(filepath, 'cwl')
        mainstr, substr_dict = NextflowTranslator.translate_workflow(wf)
    
    def test_metagoflow_conditionals_rnaprediction(self):
        filepath = 'janis_core/tests/data/cwl/MetaGOflow/workflows/conditionals/rna-prediction.cwl'
        wf = ingest(filepath, 'cwl')
        mainstr, substr_dict = NextflowTranslator.translate_workflow(wf)
    
    def test_metagoflow_subworkflows_chunking(self):
        filepath = 'janis_core/tests/data/cwl/MetaGOflow/workflows/subworkflows/chunking-subwf-IPS.cwl'
        wf = ingest(filepath, 'cwl')
        mainstr, substr_dict = NextflowTranslator.translate_workflow(wf)
    
    def test_metagoflow_subworkflows_classify(self):
        filepath = 'janis_core/tests/data/cwl/MetaGOflow/workflows/subworkflows/classify-otu-visualise.cwl'
        wf = ingest(filepath, 'cwl')
        mainstr, substr_dict = NextflowTranslator.translate_workflow(wf)
    
    def test_metagoflow_subworkflows_final_chunking(self):
        filepath = 'janis_core/tests/data/cwl/MetaGOflow/workflows/subworkflows/final_chunking.cwl'
        wf = ingest(filepath, 'cwl')
        mainstr, substr_dict = NextflowTranslator.translate_workflow(wf)
    
    def test_metagoflow_subworkflows_func_summaries(self):
        filepath = 'janis_core/tests/data/cwl/MetaGOflow/workflows/subworkflows/func_summaries.cwl'
        wf = ingest(filepath, 'cwl')
        mainstr, substr_dict = NextflowTranslator.translate_workflow(wf)
    
    def test_metagoflow_subworkflows_other_ncrnas(self):
        filepath = 'janis_core/tests/data/cwl/MetaGOflow/workflows/subworkflows/other_ncrnas.cwl'
        wf = ingest(filepath, 'cwl')
        mainstr, substr_dict = NextflowTranslator.translate_workflow(wf)
    
    def test_metagoflow_subworkflows_rna_prediction(self):
        filepath = 'janis_core/tests/data/cwl/MetaGOflow/workflows/subworkflows/rna_prediction-sub-wf.cwl'
        wf = ingest(filepath, 'cwl')
        mainstr, substr_dict = NextflowTranslator.translate_workflow(wf)
    
    def test_metagoflow_subworkflows_amplicon_its(self):
        filepath = 'janis_core/tests/data/cwl/MetaGOflow/workflows/subworkflows/amplicon/ITS-wf.cwl'
        wf = ingest(filepath, 'cwl')
        mainstr, substr_dict = NextflowTranslator.translate_workflow(wf)
    
    def test_ebi_amplicon(self):
        filepath = 'janis_core/tests/data/cwl/ebi-metagenomics/workflows/amplicon-wf--v.5-cond.cwl'
        wf = ingest(filepath, 'cwl')
        mainstr, substr_dict = NextflowTranslator.translate_workflow(wf)
    
    def test_ebi_assembly(self):
        filepath = 'janis_core/tests/data/cwl/ebi-metagenomics/workflows/assembly-wf--v.5-cond.cwl'
        wf = ingest(filepath, 'cwl')
        mainstr, substr_dict = NextflowTranslator.translate_workflow(wf)
    
    def test_ebi_raw_reads(self):
        filepath = 'janis_core/tests/data/cwl/ebi-metagenomics/workflows/raw-reads-wf--v.5-cond.cwl'
        wf = ingest(filepath, 'cwl')
        mainstr, substr_dict = NextflowTranslator.translate_workflow(wf)
    
    def test_ebi_subwf_chunking_subwf_hmmsearch(self):
        filepath = 'janis_core/tests/data/cwl/ebi-metagenomics/workflows/subworkflows/chunking-subwf-hmmsearch.cwl'
        wf = ingest(filepath, 'cwl')
        mainstr, substr_dict = NextflowTranslator.translate_workflow(wf)
    
    def test_ebi_subwf_chunking_subwf_ips(self):
        filepath = 'janis_core/tests/data/cwl/ebi-metagenomics/workflows/subworkflows/chunking-subwf-IPS.cwl'
        wf = ingest(filepath, 'cwl')
        mainstr, substr_dict = NextflowTranslator.translate_workflow(wf)
