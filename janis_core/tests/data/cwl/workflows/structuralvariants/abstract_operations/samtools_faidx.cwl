cwlVersion: v1.2
class: Operation
id: samtools_faidx
label: samtools_faidx

requirements:
  DockerRequirement:
    dockerPull: 'quay.io/biocontainers/samtools:1.5--2'

inputs:
  input_fasta:
    type: File
  input_index:
    type: 'File?'

outputs:
  fai:
    type: File
