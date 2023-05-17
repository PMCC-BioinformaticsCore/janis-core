version development

import "SamToolsMpileup_1_9_0.wdl" as S

workflow samtools_mpileup_subpipeline {
  input {
    File vcf
    File bam
    File bam_bai
    File reference
    File reference_fai
    File reference_amb
    File reference_ann
    File reference_bwt
    File reference_pac
    File reference_sa
    File reference_dict
    Boolean? samtools_mpileup_countOrphans = true
    Boolean? samtools_mpileup_noBAQ = true
    Int? samtools_mpileup_minBQ = 0
    Int? samtools_mpileup_maxDepth = 10000
  }
  call S.SamToolsMpileup as samtools_mpileup {
    input:
      countOrphans=select_first([samtools_mpileup_countOrphans, true]),
      noBAQ=select_first([samtools_mpileup_noBAQ, true]),
      maxDepth=select_first([samtools_mpileup_maxDepth, 10000]),
      positions=vcf,
      minBQ=select_first([samtools_mpileup_minBQ, 0]),
      reference=reference,
      bam=bam,
      bam_bai=bam_bai
  }
  output {
    File out = samtools_mpileup.out
  }
}