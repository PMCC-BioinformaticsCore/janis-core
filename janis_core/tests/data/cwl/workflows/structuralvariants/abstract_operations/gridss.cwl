cwlVersion: v1.2
class: Operation
id: gridss
label: gridss

doc: |
  Genomic Rearrangement IDentification Software Suite

requirements:
  DockerRequirement:
    dockerPull: 'quay.io/biocontainers/gridss:2.9.3--0'

inputs:
  input:
    type: File
    secondaryFiles:
      - .bai
  reference_genome:
    type: File
  assemblyFilename:
    type:
      - 'string?'
      - 'null'
  threads:
    type: 'int?'
  blacklist:
    type: 'File?'

outputs:
  output:
    type: File