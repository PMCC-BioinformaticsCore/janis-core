#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0
label: Filter Vardict Somatic Vcf

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: biocontainers/bcftools:v1.9-1-deb_cv1

inputs:
- id: vcf
  label: vcf
  type:
  - File
  - 'null'
  inputBinding:
    position: 1
- id: outputFilename
  label: outputFilename
  type:
  - string
  - 'null'
  default: generated.filter.vcf
  inputBinding:
    prefix: -o
    position: 3
    shellQuote: false

outputs:
- id: out
  label: out
  type: File
  outputBinding:
    glob: generated.filter.vcf
    loadContents: false
stdout: _stdout
stderr: _stderr
arguments:
- position: 0
  valueFrom: "bcftools filter -e 'STATUS=\\\"GERMLINE\\\"' -o - "
  shellQuote: false
- position: 2
  valueFrom: "| bcftools filter -i 'FILTER==\\\"PASS\\\"'"
  shellQuote: false
id: FilterVardictSomaticVcf
