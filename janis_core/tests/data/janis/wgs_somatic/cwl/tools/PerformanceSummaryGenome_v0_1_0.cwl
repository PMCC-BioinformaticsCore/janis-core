#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.0
label: Performance summary workflow (whole genome)

requirements:
- class: InlineJavascriptRequirement
- class: StepInputExpressionRequirement

inputs:
- id: bam
  type: File
  secondaryFiles:
  - .bai
- id: sample_name
  type: string
- id: genome_file
  type: File
- id: samtoolsview_doNotOutputAlignmentsWithBitsSet
  doc: |-
    Do not output alignments with any bits set in INT present in the FLAG field. INT can be specified in hex by beginning with `0x' (i.e. /^0x[0-9A-F]+/) or in octal by beginning with `0' (i.e. /^0[0-7]+/) [0].
  type: string
  default: '0x400'
- id: performancesummary_genome
  doc: |-
    calculate statistics for whole genome data.--target_flagstat must not be speicified
  type: boolean
  default: true

outputs:
- id: performanceSummaryOut
  type: File
  outputSource: performancesummary/out

steps:
- id: gatk4collectinsertsizemetrics
  label: 'GATK4: CollectInsertSizeMetrics'
  in:
  - id: bam
    source: bam
  run: Gatk4CollectInsertSizeMetrics_4_1_3_0.cwl
  out:
  - id: out
  - id: outHistogram
- id: bamflagstat
  label: 'SamTools: Flagstat'
  in:
  - id: bam
    source: bam
  run: SamToolsFlagstat_1_9_0.cwl
  out:
  - id: out
- id: samtoolsview
  label: 'SamTools: View'
  in:
  - id: doNotOutputAlignmentsWithBitsSet
    source: samtoolsview_doNotOutputAlignmentsWithBitsSet
  - id: sam
    source: bam
  run: SamToolsView_1_9_0.cwl
  out:
  - id: out
- id: rmdupbamflagstat
  label: 'SamTools: Flagstat'
  in:
  - id: bam
    source: samtoolsview/out
  run: SamToolsFlagstat_1_9_0.cwl
  out:
  - id: out
- id: bedtoolsgenomecoveragebed
  label: 'BEDTools: genomeCoverageBed'
  in:
  - id: inputBam
    source: samtoolsview/out
  - id: genome
    source: genome_file
  run: bedtoolsgenomeCoverageBed_v2_29_2.cwl
  out:
  - id: out
- id: performancesummary
  label: Performance Summary
  in:
  - id: flagstat
    source: bamflagstat/out
  - id: collectInsertSizeMetrics
    source: gatk4collectinsertsizemetrics/out
  - id: coverage
    source: bedtoolsgenomecoveragebed/out
  - id: outputPrefix
    source: sample_name
  - id: rmdupFlagstat
    source: rmdupbamflagstat/out
  - id: genome
    source: performancesummary_genome
  run: performanceSummary_0_0_7.cwl
  out:
  - id: out
id: PerformanceSummaryGenome
