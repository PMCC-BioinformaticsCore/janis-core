#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0
label: Combine Variants
doc: |2

  usage: combine_vcf.py [-h] -i I --columns COLUMNS -o O --type
                        {germline,somatic} [--regions REGIONS] [--normal NORMAL]
                        [--tumor TUMOR] [--priority PRIORITY [PRIORITY ...]]

  Extracts and combines the information from germline / somatic vcfs into one

  required arguments:
    -i I                  input vcfs, the priority of the vcfs will be based on
                          the order of the input. This parameter can be
                          specified more than once
    --columns COLUMNS     Columns to keep. This parameter can be specified more
                          than once
    -o O                  output vcf (unsorted)
    --type {germline,somatic}
                          must be either germline or somatic
    --regions REGIONS     Region file containing all the variants, used as
                          samtools mpileup
    --normal NORMAL       Sample id of germline vcf, or normal sample id of
                          somatic vcf
    --tumor TUMOR         tumor sample ID, required if inputs are somatic vcfs
    --priority PRIORITY [PRIORITY ...]
                          The priority of the callers, must match with the
                          callers in the source header

  optional arguments:
    -h, --help            show this help message and exit

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: michaelfranklin/pmacutil:0.0.8

inputs:
- id: outputFilename
  label: outputFilename
  type:
  - string
  - 'null'
  default: generated.combined.vcf
  inputBinding:
    prefix: -o
- id: vcfs
  label: vcfs
  doc: input vcfs, the priority of the vcfs will be based on the order of the input
  type:
    type: array
    inputBinding:
      prefix: -i
    items: File
  inputBinding: {}
- id: type
  label: type
  doc: germline | somatic
  type: string
  inputBinding:
    prefix: --type
- id: columns
  label: columns
  doc: Columns to keep, seperated by space output vcf (unsorted)
  type:
  - type: array
    items: string
  - 'null'
  inputBinding:
    prefix: --columns
    itemSeparator: ','
- id: normal
  label: normal
  doc: Sample id of germline vcf, or normal sample id of somatic vcf
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --normal
- id: tumor
  label: tumor
  doc: tumor sample ID, required if inputs are somatic vcfs
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --tumor
- id: priority
  label: priority
  doc: The priority of the callers, must match with the callers in the source header
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --priority

outputs:
- id: out
  label: out
  type: File
  outputBinding:
    glob: generated.combined.vcf
    loadContents: false
stdout: _stdout
stderr: _stderr

baseCommand: combine_vcf.py
arguments: []
id: combinevariants
