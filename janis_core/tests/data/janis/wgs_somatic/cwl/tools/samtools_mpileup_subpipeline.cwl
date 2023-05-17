#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.0

requirements:
- class: InlineJavascriptRequirement
- class: StepInputExpressionRequirement

inputs:
- id: vcf
  type: File
- id: bam
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
- id: samtools_mpileup_countOrphans
  doc: do not discard anomalous read pairs
  type: boolean
  default: true
- id: samtools_mpileup_noBAQ
  doc: disable BAQ (per-Base Alignment Quality)
  type: boolean
  default: true
- id: samtools_mpileup_minBQ
  doc: Minimum base quality for a base to be considered [13]
  type: int
  default: 0
- id: samtools_mpileup_maxDepth
  doc: max per-file depth; avoids excessive memory usage [8000]
  type: int
  default: 10000

outputs:
- id: out
  type: File
  outputSource: samtools_mpileup/out

steps:
- id: samtools_mpileup
  label: 'SamTools: Mpileup'
  in:
  - id: countOrphans
    source: samtools_mpileup_countOrphans
  - id: noBAQ
    source: samtools_mpileup_noBAQ
  - id: maxDepth
    source: samtools_mpileup_maxDepth
  - id: positions
    source: vcf
  - id: minBQ
    source: samtools_mpileup_minBQ
  - id: reference
    source: reference
  - id: bam
    source: bam
  run: SamToolsMpileup_1_9_0.cwl
  out:
  - id: out
id: samtools_mpileup_subpipeline
