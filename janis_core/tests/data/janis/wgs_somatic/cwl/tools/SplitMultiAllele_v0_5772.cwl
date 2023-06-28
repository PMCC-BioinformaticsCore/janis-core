#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0
label: Split Multiple Alleles

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: heuermh/vt

inputs:
- id: vcf
  label: vcf
  type: File
  inputBinding:
    position: 1
    shellQuote: false
- id: reference
  label: reference
  type: File
  secondaryFiles:
  - .fai
  - .amb
  - .ann
  - .bwt
  - .pac
  - .sa
  - ^.dict
  inputBinding:
    prefix: -r
    position: 4
    shellQuote: false
- id: outputFilename
  label: outputFilename
  type:
  - string
  - 'null'
  default: generated.norm.vcf
  inputBinding:
    prefix: -o
    position: 6
    shellQuote: false

outputs:
- id: out
  label: out
  type: File
  outputBinding:
    glob: generated.norm.vcf
    loadContents: false
stdout: _stdout
stderr: _stderr
arguments:
- position: 0
  valueFrom: 'vt decompose -s '
  shellQuote: false
- position: 2
  valueFrom: '| vt normalize -n -q - '
  shellQuote: false
id: SplitMultiAllele
