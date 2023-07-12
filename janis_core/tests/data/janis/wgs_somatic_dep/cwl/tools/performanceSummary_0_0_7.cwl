#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0
label: Performance Summary
doc: |-
  usage: performance_summary.py [-h] --flagstat FLAGSTAT
                                --collect_insert_metrics COLLECT_INSERT_METRICS
                                --coverage COVERAGE -o O
                                [--target_flagstat TARGET_FLAGSTAT]
                                [--rmdup_flagstat RMDUP_FLAGSTAT] [--genome]

  Performance summary of bam

  required arguments:
    --flagstat FLAGSTAT   output of samtools flagstat on bam
    --collect_insert_metrics COLLECT_INSERT_METRICS
                          output of CollectInsertMetrics (GATK or Picard) on bam
    --coverage COVERAGE   output of bedtools coverageBed for targeted bam;
                          bedtools genomeCoverageBed for whole genome bam
    -o O                  output summary csv name

  optional arguments:
    -h, --help            show this help message and exit
    --target_flagstat TARGET_FLAGSTAT
                          output of samtools flagstat of bam target on target
                          bed. Only specified for targeted bam
    --rmdup_flagstat RMDUP_FLAGSTAT
                          output of samtools flagstat of removed duplicates bam.
                          File to be used to extract mapping infomation if
                          specified, instead of the --flagstat file.
    --genome              calculate statistics for whole genome data.
                          --target_flagstat must not be speicified
          

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: michaelfranklin/pmacutil:0.0.7

inputs:
- id: flagstat
  label: flagstat
  doc: output of samtools flagstat on bam
  type: File
  inputBinding:
    prefix: --flagstat
- id: collectInsertSizeMetrics
  label: collectInsertSizeMetrics
  doc: output of CollectInsertMetrics (GATK or Picard) on bam
  type: File
  inputBinding:
    prefix: --collect_insert_metrics
- id: coverage
  label: coverage
  doc: |-
    output of bedtools coverageBed for targeted bam; bedtools genomeCoverageBed for whole genome bam
  type: File
  inputBinding:
    prefix: --coverage
- id: outputPrefix
  label: outputPrefix
  doc: prefix of output summary csv
  type:
  - string
  - 'null'
  default: generated.csv
  inputBinding:
    prefix: -o
- id: targetFlagstat
  label: targetFlagstat
  doc: |-
    output of samtools flagstat of bam target on target bed. Only specified for targeted bam
  type:
  - File
  - 'null'
  inputBinding:
    prefix: --target_flagstat
- id: rmdupFlagstat
  label: rmdupFlagstat
  doc: |-
    output of samtools flagstat of removed duplicates bam. File to be used to extract mapping infomation if specified, instead of the --flagstat file.
  type:
  - File
  - 'null'
  inputBinding:
    prefix: --rmdup_flagstat
- id: genome
  label: genome
  doc: |-
    calculate statistics for whole genome data.--target_flagstat must not be speicified
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --genome

outputs:
- id: out
  label: out
  type: File
  outputBinding:
    glob: $((inputs.outputPrefix + ".csv"))
    outputEval: $((inputs.outputPrefix.basename + ".csv"))
    loadContents: false
stdout: _stdout
stderr: _stderr

baseCommand: performance_summary.py
arguments: []
id: performanceSummary
