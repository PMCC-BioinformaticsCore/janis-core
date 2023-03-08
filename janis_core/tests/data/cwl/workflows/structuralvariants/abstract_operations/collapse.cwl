cwlVersion: v1.2
class: Operation
id: collapse
label: collapse

doc: |
  Merge CNVs predictions for each tool

requirements:
  DockerRequirement:
    dockerPull: 'quay.io/biocontainers/bedtools:2.26.0gx--he513fc3_4'

inputs:
  input:
    type: File

outputs:
  output:
    type: File