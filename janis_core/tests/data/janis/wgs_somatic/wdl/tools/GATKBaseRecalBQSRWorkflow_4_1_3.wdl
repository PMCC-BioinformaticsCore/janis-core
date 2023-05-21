version development

import "Gatk4BaseRecalibrator_4_1_3_0.wdl" as G
import "Gatk4ApplyBQSR_4_1_3_0.wdl" as G2

workflow GATKBaseRecalBQSRWorkflow {
  input {
    File bam
    File bam_bai
    File? intervals
    File reference
    File reference_fai
    File reference_amb
    File reference_ann
    File reference_bwt
    File reference_pac
    File reference_sa
    File reference_dict
    File snps_dbsnp
    File snps_dbsnp_tbi
    File snps_1000gp
    File snps_1000gp_tbi
    File known_indels
    File known_indels_tbi
    File mills_indels
    File mills_indels_tbi
  }
  call G.Gatk4BaseRecalibrator as base_recalibrator {
    input:
      bam=bam,
      bam_bai=bam_bai,
      knownSites=[snps_dbsnp, snps_1000gp, known_indels, mills_indels],
      knownSites_tbi=[snps_dbsnp_tbi, snps_1000gp_tbi, known_indels_tbi, mills_indels_tbi],
      reference=reference,
      reference_fai=reference_fai,
      reference_amb=reference_amb,
      reference_ann=reference_ann,
      reference_bwt=reference_bwt,
      reference_pac=reference_pac,
      reference_sa=reference_sa,
      reference_dict=reference_dict,
      intervals=intervals
  }
  call G2.Gatk4ApplyBQSR as apply_bqsr {
    input:
      bam=bam,
      bam_bai=bam_bai,
      reference=reference,
      reference_fai=reference_fai,
      reference_amb=reference_amb,
      reference_ann=reference_ann,
      reference_bwt=reference_bwt,
      reference_pac=reference_pac,
      reference_sa=reference_sa,
      reference_dict=reference_dict,
      recalFile=base_recalibrator.out,
      intervals=intervals
  }
  output {
    File out = apply_bqsr.out
    File out_bai = apply_bqsr.out_bai
  }
}