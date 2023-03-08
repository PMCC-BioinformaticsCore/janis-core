cwlVersion: v1.2
class: Operation
id: tabix
label: tabix

doc: |
  Indexing BED file

requirements:
  DockerRequirement:
    dockerPull: 'quay.io/biocontainers/tabix:0.2.6--ha92aebf_0'

inputs:
  input:
    type: File

outputs:
  output:
    type: File
    secondaryFiles: .tbi