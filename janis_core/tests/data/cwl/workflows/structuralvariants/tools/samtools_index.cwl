cwlVersion: v1.0
class: CommandLineTool
id: samtools_index
label: samtools_index

doc: |
  Indexing BAM

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - $(inputs.input)

hints:
  DockerRequirement:
    dockerPull: 'quay.io/biocontainers/samtools:1.5--2'

baseCommand: [samtools, index]

arguments:
  - position: 1
    valueFrom: -b

inputs:
  input:
    type: File
    inputBinding:
      position: 2

outputs:
  output:
    type: File
    secondaryFiles: .bai
    outputBinding:
      glob: '*.sorted*.bam'
