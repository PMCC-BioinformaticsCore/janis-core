cwlVersion: v1.0
class: CommandLineTool
id: samtools_view
label: samtools_view

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - $(inputs.input)

hints:
  DockerRequirement:
    dockerPull: 'quay.io/biocontainers/samtools:1.5--2'

baseCommand: [samtools, view]

arguments:
  - position: 5
    prefix: '-o'
    valueFrom: $(inputs.input.nameroot).filtered.bam

inputs:
  input:
    type: File
    inputBinding:
      position: 4
  min_mapping_quality:
    type: int
    inputBinding:
      position: 2
      prefix: '-q'
  bits_set:
    type: int
    inputBinding:
      position: 3
      prefix: '-F'
  threads:
    type: int?
    default: 8
    inputBinding:
      position: 1
      prefix: '--threads'

outputs:
  output:
    type: File
    outputBinding:
      glob: '*.filtered.bam'
