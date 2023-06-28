cwlVersion: v1.2
class: Operation
id: samtools_merge
label: samtools_merge

doc: |
  Merge multiple BAM files

requirements:
  DockerRequirement:
    dockerPull: 'quay.io/biocontainers/samtools:1.5--2'

inputs:
  paired:
    type: File
  unpaired_R1:
    type: File
  unpaired_R2:
    type: File
  threads:
    type: 'int?'

outputs:
  output:
    type: File