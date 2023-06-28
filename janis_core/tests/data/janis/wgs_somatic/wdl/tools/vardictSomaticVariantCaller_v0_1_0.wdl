version development

import "vardict_somatic_1_6_0.wdl" as V
import "bcftoolsAnnotate_v1_5.wdl" as B
import "bgzip_1_2_1.wdl" as B2
import "tabix_1_2_1.wdl" as T
import "SplitMultiAllele_v0_5772.wdl" as S
import "trimIUPAC_0_0_5.wdl" as T2
import "FilterVardictSomaticVcf_v1_9.wdl" as F

workflow vardictSomaticVariantCaller {
  input {
    File normal_bam
    File normal_bam_bai
    File tumor_bam
    File tumor_bam_bai
    String normal_name
    String tumor_name
    File intervals
    Float? allele_freq_threshold = 0.05
    File header_lines
    File reference
    File reference_fai
    File reference_amb
    File reference_ann
    File reference_bwt
    File reference_pac
    File reference_sa
    File reference_dict
    Boolean? vardict_chromNamesAreNumbers = true
    Boolean? vardict_vcfFormat = true
    Int? vardict_chromColumn = 1
    Int? vardict_regStartCol = 2
    Int? vardict_geneEndCol = 3
    Boolean? compressvcf_stdout = true
  }
  call V.vardict_somatic as vardict {
    input:
      tumorBam=tumor_bam,
      tumorBam_bai=tumor_bam_bai,
      normalBam=normal_bam,
      normalBam_bai=normal_bam_bai,
      intervals=intervals,
      reference=reference,
      reference_fai=reference_fai,
      tumorName=tumor_name,
      normalName=normal_name,
      alleleFreqThreshold=select_first([allele_freq_threshold, 0.05]),
      chromNamesAreNumbers=select_first([vardict_chromNamesAreNumbers, true]),
      chromColumn=select_first([vardict_chromColumn, 1]),
      geneEndCol=select_first([vardict_geneEndCol, 3]),
      regStartCol=select_first([vardict_regStartCol, 2]),
      vcfFormat=select_first([vardict_vcfFormat, true])
  }
  call B.bcftoolsAnnotate as annotate {
    input:
      vcf=vardict.out,
      headerLines=header_lines
  }
  call B2.bgzip as compressvcf {
    input:
      file=annotate.out,
      stdout=select_first([compressvcf_stdout, true])
  }
  call T.tabix as tabixvcf {
    input:
      inp=compressvcf.out
  }
  call S.SplitMultiAllele as splitnormalisevcf {
    input:
      vcf=annotate.out,
      reference=reference,
      reference_fai=reference_fai,
      reference_amb=reference_amb,
      reference_ann=reference_ann,
      reference_bwt=reference_bwt,
      reference_pac=reference_pac,
      reference_sa=reference_sa,
      reference_dict=reference_dict
  }
  call T2.trimIUPAC as trim {
    input:
      vcf=splitnormalisevcf.out
  }
  call F.FilterVardictSomaticVcf as filterpass {
    input:
      vcf=trim.out
  }
  output {
    File variants = tabixvcf.out
    File variants_tbi = tabixvcf.out_tbi
    File out = filterpass.out
  }
}