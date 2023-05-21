#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0
label: Trim IUPAC Bases

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: michaelfranklin/pmacutil:0.0.5

inputs:
- id: vcf
  label: vcf
  doc: The VCF to remove the IUPAC bases from
  type: File
  inputBinding:
    position: 0
- id: outputFilename
  label: outputFilename
  type:
  - string
  - 'null'
  default: generated.trimmed.vcf
  inputBinding:
    position: 2

outputs:
- id: out
  label: out
  type: File
  outputBinding:
    glob: generated.trimmed.vcf
    loadContents: false
stdout: _stdout
stderr: _stderr

baseCommand: trimIUPAC.py
arguments: []
id: trimIUPAC
