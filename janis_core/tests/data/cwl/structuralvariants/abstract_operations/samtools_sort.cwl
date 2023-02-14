cwlVersion: v1.2
class: Operation
id: samtools_sort
label: samtools_sort

doc: |
  Sort a bam file by read names

requirements:
  DockerRequirement:
    dockerPull: 'quay.io/biocontainers/samtools:1.5--2'

inputs:
  input:
    type: File
  threads:
    type: 'int?'

outputs:
  output:
    type: File