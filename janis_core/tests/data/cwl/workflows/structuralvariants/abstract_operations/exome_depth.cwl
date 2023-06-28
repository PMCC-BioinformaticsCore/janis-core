cwlVersion: v1.2
class: Operation
id: exome_depth
label: exome_depth

doc: |
  ExomeDepth

requirements:
  DockerRequirement:
    dockerPull: 'quay.io/biocontainers/r-exomedepth:1.1.12--r36h6786f55_0'

inputs:
  input:
    type: 'File[]'
  mapping:
    type: File
  reference_genome:
    type: File
    
outputs:
  output:
    type: File