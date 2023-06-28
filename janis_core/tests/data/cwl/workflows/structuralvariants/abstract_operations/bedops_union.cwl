cwlVersion: v1.2
class: Operation
id: bedops_union
label: bedops_union

doc: |
  Concat BED files

requirements:
  DockerRequirement:
    dockerPull: 'quay.io/biocontainers/bedops:2.4.39--h7d875b9_1'

inputs:
  input:
    type: 'File[]'

outputs:
  output:
    type: File