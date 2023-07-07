version development

import "Gatk4SplitReads_4_1_3_0.wdl" as G

workflow split_bam_subpipeline {
  input {
    File bam
    File bam_bai
    File? intervals
  }
  call G.Gatk4SplitReads as split_bam {
    input:
      bam=bam,
      bam_bai=bam_bai,
      intervals=intervals
  }
  output {
    File out = split_bam.out
    File out_bai = split_bam.out_bai
  }
}