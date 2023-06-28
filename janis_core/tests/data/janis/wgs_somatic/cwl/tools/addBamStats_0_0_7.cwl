#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0
label: Add Bam Statistics to Vcf
doc: |-
  usage: add_bam_stats.py [-h] -i I -o O --type {germline,somatic}
                          [--mpileup MPILEUP] [--normal_mpileup NORMAL_MPILEUP]
                          [--tumor_mpileup TUMOR_MPILEUP]
                          [--normal_id NORMAL_ID] [--tumor_id TUMOR_ID]

  Get stats from bam file and write to vcf

  required arguments:
    -i I                  input vcf
    -o O                  output vcf
    --type {germline,somatic}
                          must be either germline or somatic
    --mpileup MPILEUP     mpileup file extracted from bam file
    --normal_mpileup NORMAL_MPILEUP
                          mpileup file extracted from the normal sample bam,
                          required if input is somatic vcf
    --tumor_mpileup TUMOR_MPILEUP
                          mpileup file extracted from the tumor sample, required
                          if input is somatic vcf
    --normal_id NORMAL_ID
                          Normal sample id, required if input is somatic vcf
    --tumor_id TUMOR_ID   Tumor sample id, required if input is somatic vcf

  optional arguments:
    -h, --help            show this help message and exit
          

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: michaelfranklin/pmacutil:0.0.7

inputs:
- id: mpileup
  label: mpileup
  doc: mpileup file extracted from bam file
  type:
  - File
  - 'null'
  inputBinding:
    prefix: --mpileup
- id: normalMpileup
  label: normalMpileup
  doc: |-
    mpileup file extracted from the normal sample bam, required if input is somatic vcf
  type:
  - File
  - 'null'
  inputBinding:
    prefix: --normal_mpileup
- id: tumorMpileup
  label: tumorMpileup
  doc: |-
    mpileup file extracted from the tumor sample bam, required if input is somatic vcf
  type:
  - File
  - 'null'
  inputBinding:
    prefix: --tumor_mpileup
- id: normalID
  label: normalID
  doc: normal sample id, required if input is somatic vcf
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --normal_id
- id: tumorID
  label: tumorID
  doc: tumor sample id, required if input is somatic vcf
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --tumor_id
- id: inputVcf
  label: inputVcf
  doc: input vcf
  type: File
  inputBinding:
    prefix: -i
- id: outputFilename
  label: outputFilename
  doc: output vcf name
  type:
  - string
  - 'null'
  default: generated.addbamstats.vcf
  inputBinding:
    prefix: -o
- id: type
  label: type
  doc: must be either germline or somatic
  type: string
  inputBinding:
    prefix: --type

outputs:
- id: out
  label: out
  type: File
  outputBinding:
    glob: generated.addbamstats.vcf
    loadContents: false
stdout: _stdout
stderr: _stderr

baseCommand: add_bam_stats.py
arguments: []
id: addBamStats
