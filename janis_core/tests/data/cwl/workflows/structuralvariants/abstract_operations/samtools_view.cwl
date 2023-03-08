cwlVersion: v1.2
class: Operation
id: samtools_view
label: samtools_view

requirements:
  DockerRequirement:
    dockerPull: 'quay.io/biocontainers/samtools:1.5--2'

inputs:
  input:
    type: File
  min_mapping_quality:
    type: int
  bits_set:
    type: int
  threads:
    type: 'int?'

outputs:
  output:
    type: File