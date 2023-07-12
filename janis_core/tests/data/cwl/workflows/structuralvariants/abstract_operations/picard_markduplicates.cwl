cwlVersion: v1.2
class: Operation
id: picard_markduplicates
label: picard_markduplicates

doc: |
  Removal of duplicates from aligned reads

requirements:
  DockerRequirement:
    dockerPull: 'quay.io/biocontainers/picard:2.10.6--py35_0'

inputs:
  input:
    type: File

outputs:
  alignments:
    type: File
  metrics:
    type: File