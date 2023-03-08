cwlVersion: v1.0
class: CommandLineTool
id: samtools_faidx
label: samtools_faidx

requirements:
  InlineJavascriptRequirement: {}
  ShellCommandRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - $(inputs.input_fasta)
      - $(inputs.input_index)

hints:
  DockerRequirement:
    dockerPull: 'quay.io/biocontainers/samtools:1.5--2'

baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      $(inputs.input_index ? 'echo samtools faidx' : 'samtools faidx')

inputs:
  input_fasta:
    type: File
    inputBinding:
      position: 1
  input_index:
    type: File?

outputs:
  fai:
    type: File
    outputBinding:
      glob: '*.fai'
