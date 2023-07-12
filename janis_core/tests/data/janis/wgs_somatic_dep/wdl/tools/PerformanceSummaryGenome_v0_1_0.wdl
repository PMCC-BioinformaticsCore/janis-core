version development

import "Gatk4CollectInsertSizeMetrics_4_1_3_0.wdl" as G
import "SamToolsFlagstat_1_9_0.wdl" as S
import "SamToolsView_1_9_0.wdl" as S2
import "bedtoolsgenomeCoverageBed_v2_29_2.wdl" as B
import "performanceSummary_0_0_7.wdl" as P

workflow PerformanceSummaryGenome {
  input {
    File bam
    File bam_bai
    String sample_name
    File genome_file
    String? samtoolsview_doNotOutputAlignmentsWithBitsSet = "0x400"
    Boolean? performancesummary_genome = true
  }
  call G.Gatk4CollectInsertSizeMetrics as gatk4collectinsertsizemetrics {
    input:
      bam=bam,
      bam_bai=bam_bai
  }
  call S.SamToolsFlagstat as bamflagstat {
    input:
      bam=bam
  }
  call S2.SamToolsView as samtoolsview {
    input:
      doNotOutputAlignmentsWithBitsSet=select_first([samtoolsview_doNotOutputAlignmentsWithBitsSet, "0x400"]),
      sam=bam
  }
  call S.SamToolsFlagstat as rmdupbamflagstat {
    input:
      bam=samtoolsview.out
  }
  call B.bedtoolsgenomeCoverageBed as bedtoolsgenomecoveragebed {
    input:
      inputBam=samtoolsview.out,
      genome=genome_file
  }
  call P.performanceSummary as performancesummary {
    input:
      flagstat=bamflagstat.out,
      collectInsertSizeMetrics=gatk4collectinsertsizemetrics.out,
      coverage=bedtoolsgenomecoveragebed.out,
      outputPrefix=sample_name,
      rmdupFlagstat=rmdupbamflagstat.out,
      genome=select_first([performancesummary_genome, true])
  }
  output {
    File performanceSummaryOut = performancesummary.out
  }
}