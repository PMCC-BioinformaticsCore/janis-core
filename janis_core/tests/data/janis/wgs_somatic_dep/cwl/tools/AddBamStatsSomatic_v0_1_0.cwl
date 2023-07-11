#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.0
label: Annotate Bam Stats to Somatic Vcf Workflow

requirements:
- class: InlineJavascriptRequirement
- class: StepInputExpressionRequirement
- class: SubworkflowFeatureRequirement

inputs:
- id: normal_id
  type: string
- id: tumor_id
  type: string
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
- id: vcf
  type: File
- id: addbamstats_type
  doc: must be either germline or somatic
  type: string
  default: somatic

outputs:
- id: out
  type: File
  outputSource: addbamstats/out

steps:
- id: tumor
  in:
  - id: vcf
    source: vcf
  - id: bam
    source: tumor_bam
  - id: reference
    source: reference
  run: samtools_mpileup_subpipeline.cwl
  out:
  - id: out
- id: normal
  in:
  - id: vcf
    source: vcf
  - id: bam
    source: normal_bam
  - id: reference
    source: reference
  run: samtools_mpileup_subpipeline.cwl
  out:
  - id: out
- id: addbamstats
  label: Add Bam Statistics to Vcf
  in:
  - id: normalMpileup
    source: normal/out
  - id: tumorMpileup
    source: tumor/out
  - id: normalID
    source: normal_id
  - id: tumorID
    source: tumor_id
  - id: inputVcf
    source: vcf
  - id: type
    source: addbamstats_type
  run: addBamStats_0_0_7.cwl
  out:
  - id: out
id: AddBamStatsSomatic
