cwlVersion: v1.0
class: CommandLineTool
id: bwa_index
label: bwa_index

requirements:
  InlineJavascriptRequirement: {}
  ShellCommandRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - $(inputs.input_fasta)
      - $(inputs.input_amb)
      - $(inputs.input_ann)
      - $(inputs.input_bwt)
      - $(inputs.input_pac)
      - $(inputs.input_sa)

hints:
  DockerRequirement:
    dockerPull: 'quay.io/biocontainers/bwa:0.7.17--h84994c4_5'

baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      $(inputs.input_amb && inputs.input_ann && inputs.input_bwt && inputs.input_pac && inputs.input_sa ? 'echo bwa' : inputs.generate_bwa_indexes ? 'bwa index' : 'echo bwa')

inputs:
  generate_bwa_indexes:
    type: boolean?
  input_fasta:
    type: File
    inputBinding:
      position: 2
  input_amb:
    type: File?
  input_ann:
    type: File?
  input_bwt:
    type: File?
  input_pac:
    type: File?
  input_sa:
    type: File?
  algoType:
    type: string?
    inputBinding:
      position: 1
      prefix: '-a'

outputs:
  amb:
    type: File?
    outputBinding:
      glob: '*.amb'
  ann:
    type: File?
    outputBinding:
      glob: '*.ann'
  bwt:
    type: File?
    outputBinding:
      glob: '*.bwt'
  pac:
    type: File?
    outputBinding:
      glob: '*.pac'
  sa:
    type: File?
    outputBinding:
      glob: '*.sa'
