#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.0
label: GATK Base Recalibration on Bam

requirements:
- class: InlineJavascriptRequirement
- class: StepInputExpressionRequirement
- class: MultipleInputFeatureRequirement

inputs:
- id: bam
  type: File
  secondaryFiles:
  - .bai
- id: intervals
  doc: |-
    This optional interval supports processing by regions. If this input resolves to null, then GATK will process the whole genome per each tool's spec
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
- id: snps_dbsnp
  type: File
  secondaryFiles:
  - .tbi
- id: snps_1000gp
  type: File
  secondaryFiles:
  - .tbi
- id: known_indels
  type: File
  secondaryFiles:
  - .tbi
- id: mills_indels
  type: File
  secondaryFiles:
  - .tbi

outputs:
- id: out
  type: File
  secondaryFiles:
  - .bai
  outputSource: apply_bqsr/out

steps:
- id: base_recalibrator
  label: 'GATK4: Base Recalibrator'
  in:
  - id: bam
    source: bam
  - id: knownSites
    source:
    - snps_dbsnp
    - snps_1000gp
    - known_indels
    - mills_indels
  - id: reference
    source: reference
  - id: intervals
    source: intervals
  run: Gatk4BaseRecalibrator_4_1_3_0.cwl
  out:
  - id: out
- id: apply_bqsr
  label: 'GATK4: Apply base quality score recalibration'
  in:
  - id: bam
    source: bam
  - id: reference
    source: reference
  - id: recalFile
    source: base_recalibrator/out
  - id: intervals
    source: intervals
  run: Gatk4ApplyBQSR_4_1_3_0.cwl
  out:
  - id: out
id: GATKBaseRecalBQSRWorkflow
