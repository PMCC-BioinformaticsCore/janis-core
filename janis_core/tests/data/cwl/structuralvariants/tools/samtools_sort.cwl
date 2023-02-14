cwlVersion: v1.0
class: CommandLineTool
id: samtools_sort
label: samtools_sort

doc: |
  Sort a bam file by read names

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - $(inputs.input)

hints:
  DockerRequirement:
    dockerPull: 'quay.io/biocontainers/samtools:1.5--2'

baseCommand: [samtools, sort]

arguments:
  - position: 3
    prefix: '-o'
    valueFrom: $(inputs.input.nameroot).sorted.bam

inputs:
  input:
    type: File
    inputBinding:
      position: 2
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
      glob: '*.sorted.bam'
