cwlVersion: v1.0
class: CommandLineTool
id: samtools_view_sam2bam
label: samtools_view_sam2bam

doc: |
  Convert SAM to BAM

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
  - position: 3
    prefix: '-o'
    valueFrom: $(inputs.input.nameroot).bam

inputs:
  input:
    type: File
    inputBinding:
      position: 2
      prefix: '-bS'
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
      glob: '*.bam'
