cwlVersion: v1.2
class: Operation
id: samtools_index
label: samtools_index

doc: |
  Indexing BAM

requirements:
  DockerRequirement:
    dockerPull: 'quay.io/biocontainers/samtools:1.5--2'

inputs:
  input:
    type: File

outputs:
  output:
    type: File
    secondaryFiles: .bai