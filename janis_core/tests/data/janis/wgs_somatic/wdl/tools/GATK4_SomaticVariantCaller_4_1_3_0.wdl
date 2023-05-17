version development

import "split_bam_subpipeline.wdl" as S
import "Gatk4Mutect2_4_1_3_0.wdl" as G
import "Gatk4LearnReadOrientationModel_4_1_4_0.wdl" as G2
import "Gatk4GetPileupSummaries_4_1_6_0.wdl" as G3
import "Gatk4CalculateContamination_4_1_4_0.wdl" as G4
import "Gatk4FilterMutectCalls_4_1_3_0.wdl" as G5
import "UncompressArchive_v1_0_0.wdl" as U
import "SplitMultiAllele_v0_5772.wdl" as S2
import "VcfTools_0_1_16.wdl" as V

workflow GATK4_SomaticVariantCaller {
  input {
    File normal_bam
    File normal_bam_bai
    File tumor_bam
    File tumor_bam_bai
    String? normal_name
    File? intervals
    File reference
    File reference_fai
    File reference_amb
    File reference_ann
    File reference_bwt
    File reference_pac
    File reference_sa
    File reference_dict
    File gnomad
    File gnomad_tbi
    File? panel_of_normals
    File? panel_of_normals_tbi
    Boolean? filterpass_removeFileteredAll = true
    Boolean? filterpass_recode = true
    Boolean? filterpass_recodeINFOAll = true
  }
  call S.split_bam_subpipeline as normal_split_bam {
    input:
      bam=normal_bam,
      bam_bai=normal_bam_bai,
      intervals=intervals
  }
  call S.split_bam_subpipeline as tumor_split_bam {
    input:
      bam=tumor_bam,
      bam_bai=tumor_bam_bai,
      intervals=intervals
  }
  call G.Gatk4Mutect2 as mutect2 {
    input:
      tumorBams=[tumor_split_bam.out],
      tumorBams_bai=[tumor_split_bam.out_bai],
      normalBams=[normal_split_bam.out],
      normalBams_bai=[normal_split_bam.out_bai],
      normalSample=normal_name,
      reference=reference,
      reference_fai=reference_fai,
      reference_amb=reference_amb,
      reference_ann=reference_ann,
      reference_bwt=reference_bwt,
      reference_pac=reference_pac,
      reference_sa=reference_sa,
      reference_dict=reference_dict,
      germlineResource=gnomad,
      germlineResource_tbi=gnomad_tbi,
      intervals=intervals,
      panelOfNormals=panel_of_normals,
      panelOfNormals_tbi=panel_of_normals_tbi
  }
  call G2.Gatk4LearnReadOrientationModel as learnorientationmodel {
    input:
      f1r2CountsFiles=[mutect2.f1f2r_out]
  }
  call G3.Gatk4GetPileupSummaries as getpileupsummaries {
    input:
      bam=[tumor_split_bam.out],
      bam_bai=[tumor_split_bam.out_bai],
      sites=gnomad,
      sites_tbi=gnomad_tbi,
      intervals=intervals
  }
  call G4.Gatk4CalculateContamination as calculatecontamination {
    input:
      pileupTable=getpileupsummaries.out
  }
  call G5.Gatk4FilterMutectCalls as filtermutect2calls {
    input:
      contaminationTable=calculatecontamination.contOut,
      segmentationFile=calculatecontamination.segOut,
      statsFile=mutect2.stats,
      readOrientationModel=learnorientationmodel.out,
      vcf=mutect2.out,
      vcf_tbi=mutect2.out_tbi,
      reference=reference,
      reference_fai=reference_fai,
      reference_amb=reference_amb,
      reference_ann=reference_ann,
      reference_bwt=reference_bwt,
      reference_pac=reference_pac,
      reference_sa=reference_sa,
      reference_dict=reference_dict
  }
  call U.UncompressArchive as uncompressvcf {
    input:
      file=filtermutect2calls.out
  }
  call S2.SplitMultiAllele as splitnormalisevcf {
    input:
      vcf=uncompressvcf.out,
      reference=reference,
      reference_fai=reference_fai,
      reference_amb=reference_amb,
      reference_ann=reference_ann,
      reference_bwt=reference_bwt,
      reference_pac=reference_pac,
      reference_sa=reference_sa,
      reference_dict=reference_dict
  }
  call V.VcfTools as filterpass {
    input:
      vcf=splitnormalisevcf.out,
      removeFileteredAll=select_first([filterpass_removeFileteredAll, true]),
      recode=select_first([filterpass_recode, true]),
      recodeINFOAll=select_first([filterpass_recodeINFOAll, true])
  }
  output {
    File variants = filtermutect2calls.out
    File variants_tbi = filtermutect2calls.out_tbi
    File? out_bam = mutect2.bam
    File? out_bam_bai = mutect2.bam_bai
    File out = filterpass.out
  }
}