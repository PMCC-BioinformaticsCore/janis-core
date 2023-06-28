version development

import "manta_1_5_0.wdl" as M
import "strelka_somatic_2_9_10.wdl" as S
import "ConcatStrelkaSomaticVcf_0_1_16.wdl" as C
import "bcftoolssort_v1_9.wdl" as B
import "SplitMultiAllele_v0_5772.wdl" as S2
import "extractStrelkaSomaticADDP_0_1_1.wdl" as E
import "VcfTools_0_1_16.wdl" as V

workflow strelkaSomaticVariantCaller {
  input {
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
    File? intervals
    File? intervals_tbi
    Boolean? is_exome
    Boolean? filterpass_removeFileteredAll = true
    Boolean? filterpass_recode = true
    Boolean? filterpass_recodeINFOAll = true
  }
  call M.manta as manta {
    input:
      bam=normal_bam,
      bam_bai=normal_bam_bai,
      reference=reference,
      reference_fai=reference_fai,
      tumorBam=tumor_bam,
      tumorBam_bai=tumor_bam_bai,
      exome=is_exome,
      callRegions=intervals,
      callRegions_tbi=intervals_tbi
  }
  call S.strelka_somatic as strelka {
    input:
      normalBam=normal_bam,
      normalBam_bai=normal_bam_bai,
      tumorBam=tumor_bam,
      tumorBam_bai=tumor_bam_bai,
      reference=reference,
      reference_fai=reference_fai,
      indelCandidates=[manta.candidateSmallIndels],
      indelCandidates_tbi=[manta.candidateSmallIndels_tbi],
      exome=is_exome,
      callRegions=intervals,
      callRegions_tbi=intervals_tbi
  }
  call C.ConcatStrelkaSomaticVcf as concatvcf {
    input:
      headerVcfs=[strelka.snvs, strelka.indels],
      headerVcfs_tbi=[strelka.snvs_tbi, strelka.indels_tbi],
      contentVcfs=[strelka.snvs, strelka.indels],
      contentVcfs_tbi=[strelka.snvs_tbi, strelka.indels_tbi]
  }
  call B.bcftoolssort as sortvcf {
    input:
      vcf=concatvcf.out
  }
  call S2.SplitMultiAllele as splitnormalisevcf {
    input:
      vcf=sortvcf.out,
      reference=reference,
      reference_fai=reference_fai,
      reference_amb=reference_amb,
      reference_ann=reference_ann,
      reference_bwt=reference_bwt,
      reference_pac=reference_pac,
      reference_sa=reference_sa,
      reference_dict=reference_dict
  }
  call E.extractStrelkaSomaticADDP as extractaddp {
    input:
      vcf=splitnormalisevcf.out
  }
  call V.VcfTools as filterpass {
    input:
      vcf=extractaddp.out,
      removeFileteredAll=select_first([filterpass_removeFileteredAll, true]),
      recode=select_first([filterpass_recode, true]),
      recodeINFOAll=select_first([filterpass_recodeINFOAll, true])
  }
  output {
    File sv = manta.diploidSV
    File sv_tbi = manta.diploidSV_tbi
    File variants = sortvcf.out
    File out = filterpass.out
  }
}