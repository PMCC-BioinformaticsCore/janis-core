#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0
label: Concat Strelka Somatic Vcf

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: biocontainers/vcftools:v0.1.16-1-deb_cv1

inputs:
- id: headerVcfs
  label: headerVcfs
  type:
    type: array
    items: File
  inputBinding:
    position: 1
- id: contentVcfs
  label: contentVcfs
  type:
    type: array
    items: File
  inputBinding:
    position: 4
- id: outputFilename
  label: outputFilename
  type:
  - string
  - 'null'
  default: generated.strelka.vcf
  inputBinding:
    prefix: '>'
    position: 6
    shellQuote: false

outputs:
- id: out
  label: out
  type: File
  outputBinding:
    glob: generated.strelka.vcf
    loadContents: false
stdout: _stdout
stderr: _stderr
arguments:
- position: 0
  valueFrom: vcf-merge
  shellQuote: false
- position: 2
  valueFrom: "| grep '^##' > header.vcf;"
  shellQuote: false
- position: 3
  valueFrom: vcf-concat
  shellQuote: false
- position: 5
  valueFrom: "| grep -v '^##' > content.vcf; cat header.vcf content.vcf"
  shellQuote: false
id: ConcatStrelkaSomaticVcf
