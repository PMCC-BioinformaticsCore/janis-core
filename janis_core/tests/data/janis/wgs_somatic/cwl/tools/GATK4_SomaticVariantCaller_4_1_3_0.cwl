#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.0
label: GATK4 Somatic Variant Caller
doc: |-
  This is a VariantCaller based on the GATK Best Practice pipelines. It uses the GATK4 toolkit, specifically 4.0.12.0. Takes GATK Base Recalibrated Bam as input

          It has the following steps:

          1. Mutect2
          2. LearnOrientationModel
          3. GetPileUpSummaries
          4. CalculateContamination
          5. FilterMutectCalls
          6. Split and normliase vcf
          7. Filter PASS variants

requirements:
- class: InlineJavascriptRequirement
- class: StepInputExpressionRequirement
- class: SubworkflowFeatureRequirement
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
- id: normal_name
  type:
  - string
  - 'null'
- id: intervals
  doc: |-
    This optional intervals file supports processing by regions. If this file resolves to null, then GATK will process the whole genome per each tool's spec
  type:
  - File
  - 'null'
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
- id: gnomad
  type: File
  secondaryFiles:
  - .tbi
- id: panel_of_normals
  type:
  - File
  - 'null'
  secondaryFiles:
  - .tbi
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
- id: variants
  type: File
  secondaryFiles:
  - .tbi
  outputSource: filtermutect2calls/out
- id: out_bam
  type:
  - File
  - 'null'
  secondaryFiles:
  - .bai
  outputSource: mutect2/bam
- id: out
  type: File
  outputSource: filterpass/out

steps:
- id: normal_split_bam
  in:
  - id: bam
    source: normal_bam
  - id: intervals
    source: intervals
  run: split_bam_subpipeline.cwl
  out:
  - id: out
- id: tumor_split_bam
  in:
  - id: bam
    source: tumor_bam
  - id: intervals
    source: intervals
  run: split_bam_subpipeline.cwl
  out:
  - id: out
- id: mutect2
  label: GatkMutect2
  in:
  - id: tumorBams
    source:
    - tumor_split_bam/out
    linkMerge: merge_nested
  - id: normalBams
    source:
    - normal_split_bam/out
    linkMerge: merge_nested
  - id: normalSample
    source: normal_name
  - id: reference
    source: reference
  - id: germlineResource
    source: gnomad
  - id: intervals
    source: intervals
  - id: panelOfNormals
    source: panel_of_normals
  run: Gatk4Mutect2_4_1_3_0.cwl
  out:
  - id: out
  - id: stats
  - id: f1f2r_out
  - id: bam
- id: learnorientationmodel
  label: 'GATK4: LearnReadOrientationModel'
  in:
  - id: f1r2CountsFiles
    source:
    - mutect2/f1f2r_out
    linkMerge: merge_nested
  run: Gatk4LearnReadOrientationModel_4_1_4_0.cwl
  out:
  - id: out
- id: getpileupsummaries
  label: 'GATK4: GetPileupSummaries'
  in:
  - id: bam
    source:
    - tumor_split_bam/out
    linkMerge: merge_nested
  - id: sites
    source: gnomad
  - id: intervals
    source: intervals
  run: Gatk4GetPileupSummaries_4_1_6_0.cwl
  out:
  - id: out
- id: calculatecontamination
  label: 'GATK4: CalculateContamination'
  in:
  - id: pileupTable
    source: getpileupsummaries/out
  run: Gatk4CalculateContamination_4_1_4_0.cwl
  out:
  - id: contOut
  - id: segOut
- id: filtermutect2calls
  label: 'GATK4: GetFilterMutectCalls'
  in:
  - id: contaminationTable
    source: calculatecontamination/contOut
  - id: segmentationFile
    source: calculatecontamination/segOut
  - id: statsFile
    source: mutect2/stats
  - id: readOrientationModel
    source: learnorientationmodel/out
  - id: vcf
    source: mutect2/out
  - id: reference
    source: reference
  run: Gatk4FilterMutectCalls_4_1_3_0.cwl
  out:
  - id: out
- id: uncompressvcf
  label: UncompressArchive
  in:
  - id: file
    source: filtermutect2calls/out
  run: UncompressArchive_v1_0_0.cwl
  out:
  - id: out
- id: splitnormalisevcf
  label: Split Multiple Alleles
  in:
  - id: vcf
    source: uncompressvcf/out
  - id: reference
    source: reference
  run: SplitMultiAllele_v0_5772.cwl
  out:
  - id: out
- id: filterpass
  label: VcfTools
  in:
  - id: vcf
    source: splitnormalisevcf/out
  - id: removeFileteredAll
    source: filterpass_removeFileteredAll
  - id: recode
    source: filterpass_recode
  - id: recodeINFOAll
    source: filterpass_recodeINFOAll
  run: VcfTools_0_1_16.cwl
  out:
  - id: out
id: GATK4_SomaticVariantCaller
