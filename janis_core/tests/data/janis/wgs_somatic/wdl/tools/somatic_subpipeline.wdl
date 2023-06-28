version development

import "fastqc_v0_11_8.wdl" as F
import "ParseFastqcAdaptors_v0_1_0.wdl" as P
import "BwaAligner_1_0_0.wdl" as B
import "mergeAndMarkBams_4_1_3.wdl" as M
import "GenerateGenomeFileForBedtoolsCoverage_v0_1_0.wdl" as G
import "PerformanceSummaryGenome_v0_1_0.wdl" as P2

workflow somatic_subpipeline {
  input {
    Array[Array[File]] reads
    String sample_name
    File reference
    File reference_fai
    File reference_amb
    File reference_ann
    File reference_bwt
    File reference_pac
    File reference_sa
    File reference_dict
    File? cutadapt_adapters
    Array[File] gatk_intervals
    File snps_dbsnp
    File snps_dbsnp_tbi
    File snps_1000gp
    File snps_1000gp_tbi
    File known_indels
    File known_indels_tbi
    File mills_indels
    File mills_indels_tbi
    String? align_and_sort_sortsam_tmpDir
  }
  scatter (r in reads) {
     call F.fastqc as fastqc {
      input:
        reads=r
    }
  }
  scatter (f in fastqc.datafile) {
     call P.ParseFastqcAdaptors as getfastqc_adapters {
      input:
        fastqc_datafiles=f,
        cutadapt_adaptors_lookup=cutadapt_adapters
    }
  }
  scatter (Q in zip(reads, zip(getfastqc_adapters.adaptor_sequences, getfastqc_adapters.adaptor_sequences))) {
     call B.BwaAligner as align_and_sort {
      input:
        sample_name=sample_name,
        reference=reference,
        reference_fai=reference_fai,
        reference_amb=reference_amb,
        reference_ann=reference_ann,
        reference_bwt=reference_bwt,
        reference_pac=reference_pac,
        reference_sa=reference_sa,
        reference_dict=reference_dict,
        fastq=Q.left,
        cutadapt_adapter=Q.right.right,
        cutadapt_removeMiddle3Adapter=Q.right.right,
        sortsam_tmpDir=align_and_sort_sortsam_tmpDir
    }
  }
  call M.mergeAndMarkBams as merge_and_mark {
    input:
      bams=align_and_sort.out,
      bams_bai=align_and_sort.out_bai,
      sampleName=sample_name
  }
  call G.GenerateGenomeFileForBedtoolsCoverage as calculate_performancesummary_genomefile {
    input:
      reference=reference,
      reference_dict=reference_dict
  }
  call P2.PerformanceSummaryGenome as performance_summary {
    input:
      bam=merge_and_mark.out,
      bam_bai=merge_and_mark.out_bai,
      sample_name=sample_name,
      genome_file=calculate_performancesummary_genomefile.out
  }
  output {
    File out_bam = merge_and_mark.out
    File out_bam_bai = merge_and_mark.out_bai
    Array[Array[File]] out_fastqc_reports = fastqc.out
    File out_performance_summary = performance_summary.performanceSummaryOut
  }
}