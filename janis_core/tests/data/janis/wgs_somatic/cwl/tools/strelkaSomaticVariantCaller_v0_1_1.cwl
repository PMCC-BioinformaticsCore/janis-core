#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.0
label: Strelka Somatic Variant Caller

requirements:
- class: InlineJavascriptRequirement
- class: StepInputExpressionRequirement
- class: MultipleInputFeatureRequirement

inputs:
- id: normal_bam
  type: File
  secondaryFiles:
  - .bai
- id: tumor_bam
  type: File
  secondaryFiles:
  - .bai
- id: reference
  type: File
  secondaryFiles:
  - .fai
  - .amb
  - .ann
  - .bwt
  - .pac
  - .sa
  - ^.dict
- id: intervals
  type:
  - File
  - 'null'
  secondaryFiles:
  - .tbi
- id: is_exome
  type:
  - boolean
  - 'null'
- id: filterpass_removeFileteredAll
  doc: Removes all sites with a FILTER flag other than PASS.
  type: boolean
  default: true
- id: filterpass_recode
  doc: ''
  type: boolean
  default: true
- id: filterpass_recodeINFOAll
  doc: |-
    These options can be used with the above recode options to define an INFO key name to keep in the output  file.  This  option can be used multiple times to keep more of the INFO fields. The second option is used to keep all INFO values in the original file.
  type: boolean
  default: true

outputs:
- id: sv
  type: File
  secondaryFiles:
  - .tbi
  outputSource: manta/diploidSV
- id: variants
  type: File
  outputSource: sortvcf/out
- id: out
  type: File
  outputSource: filterpass/out

steps:
- id: manta
  label: Manta
  in:
  - id: bam
    source: normal_bam
  - id: reference
    source: reference
  - id: tumorBam
    source: tumor_bam
  - id: exome
    source: is_exome
  - id: callRegions
    source: intervals
  run: manta_1_5_0.cwl
  out:
  - id: python
  - id: pickle
  - id: candidateSV
  - id: candidateSmallIndels
  - id: diploidSV
  - id: alignmentStatsSummary
  - id: svCandidateGenerationStats
  - id: svLocusGraphStats
  - id: somaticSVs
- id: strelka
  label: Strelka (Somatic)
  in:
  - id: normalBam
    source: normal_bam
  - id: tumorBam
    source: tumor_bam
  - id: reference
    source: reference
  - id: indelCandidates
    source:
    - manta/candidateSmallIndels
    linkMerge: merge_nested
  - id: exome
    source: is_exome
  - id: callRegions
    source: intervals
  run: strelka_somatic_2_9_10.cwl
  out:
  - id: configPickle
  - id: script
  - id: stats
  - id: indels
  - id: snvs
- id: concatvcf
  label: Concat Strelka Somatic Vcf
  in:
  - id: headerVcfs
    source:
    - strelka/snvs
    - strelka/indels
  - id: contentVcfs
    source:
    - strelka/snvs
    - strelka/indels
  run: ConcatStrelkaSomaticVcf_0_1_16.cwl
  out:
  - id: out
- id: sortvcf
  label: 'BCFTools: Sort'
  in:
  - id: vcf
    source: concatvcf/out
  run: bcftoolssort_v1_9.cwl
  out:
  - id: out
- id: splitnormalisevcf
  label: Split Multiple Alleles
  in:
  - id: vcf
    source: sortvcf/out
  - id: reference
    source: reference
  run: SplitMultiAllele_v0_5772.cwl
  out:
  - id: out
- id: extractaddp
  label: Extract Strelka Somatic AD DP
  in:
  - id: vcf
    source: splitnormalisevcf/out
  run: extractStrelkaSomaticADDP_0_1_1.cwl
  out:
  - id: out
- id: filterpass
  label: VcfTools
  in:
  - id: vcf
    source: extractaddp/out
  - id: removeFileteredAll
    source: filterpass_removeFileteredAll
  - id: recode
    source: filterpass_recode
  - id: recodeINFOAll
    source: filterpass_recodeINFOAll
  run: VcfTools_0_1_16.cwl
  out:
  - id: out
id: strelkaSomaticVariantCaller
