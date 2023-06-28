#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.0

requirements:
- class: InlineJavascriptRequirement
- class: StepInputExpressionRequirement

inputs:
- id: bam
  type: File
  secondaryFiles:
  - .bai
- id: intervals
  type:
  - File
  - 'null'

outputs:
- id: out
  type: File
  secondaryFiles:
  - .bai
  outputSource: split_bam/out

steps:
- id: split_bam
  label: 'GATK4: SplitReads'
  in:
  - id: bam
    source: bam
  - id: intervals
    source: intervals
  run: Gatk4SplitReads_4_1_3_0.cwl
  out:
  - id: out
id: split_bam_subpipeline
