version development

import "samtools_mpileup_subpipeline.wdl" as S
import "addBamStats_0_0_7.wdl" as A

workflow AddBamStatsSomatic {
  input {
    String normal_id
    String tumor_id
    File normal_bam
    File normal_bam_bai
    File tumor_bam
    File tumor_bam_bai
    File reference
    File reference_fai
    File reference_amb
    File reference_ann
    File reference_bwt
    File reference_pac
    File reference_sa
    File reference_dict
    File vcf
    String? addbamstats_type = "somatic"
  }
  call S.samtools_mpileup_subpipeline as tumor {
    input:
      vcf=vcf,
      bam=tumor_bam,
      bam_bai=tumor_bam_bai,
      reference=reference,
      reference_fai=reference_fai,
      reference_amb=reference_amb,
      reference_ann=reference_ann,
      reference_bwt=reference_bwt,
      reference_pac=reference_pac,
      reference_sa=reference_sa,
      reference_dict=reference_dict
  }
  call S.samtools_mpileup_subpipeline as normal {
    input:
      vcf=vcf,
      bam=normal_bam,
      bam_bai=normal_bam_bai,
      reference=reference,
      reference_fai=reference_fai,
      reference_amb=reference_amb,
      reference_ann=reference_ann,
      reference_bwt=reference_bwt,
      reference_pac=reference_pac,
      reference_sa=reference_sa,
      reference_dict=reference_dict
  }
  call A.addBamStats as addbamstats {
    input:
      normalMpileup=normal.out,
      tumorMpileup=tumor.out,
      normalID=normal_id,
      tumorID=tumor_id,
      inputVcf=vcf,
      type=select_first([addbamstats_type, "somatic"])
  }
  output {
    File out = addbamstats.out
  }
}