cwlVersion: v1.2
class: Operation
id: batch_parser
label: batch_parser

doc: |
  Group samples by batch

requirements:
  DockerRequirement:
    dockerPull: 'r-base:latest'

inputs:
  input:
    type: 'File[]'
  samples:
    type: File

outputs:
  output:
    type: 'File[]'