cwlVersion: v1.0
class: CommandLineTool
id: samtools_merge
label: samtools_merge

doc: |
  Merge multiple BAM files

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - $(inputs.paired)
      - $(inputs.unpaired_R1)
      - $(inputs.unpaired_R2)

hints:
  DockerRequirement:
    dockerPull: 'quay.io/biocontainers/samtools:1.5--2'

baseCommand: [samtools, merge, -f]

arguments:
  - position: 2
    valueFrom: $(inputs.paired.nameroot.split('.')[0]).sorted.bam

inputs:
  paired:
    type: File
    inputBinding:
      position: 3
  unpaired_R1:
    type: File
    inputBinding:
      position: 4
  unpaired_R2:
    type: File
    inputBinding:
      position: 5
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
      glob: $(inputs.paired.nameroot.split('.')[0]).sorted.bam
