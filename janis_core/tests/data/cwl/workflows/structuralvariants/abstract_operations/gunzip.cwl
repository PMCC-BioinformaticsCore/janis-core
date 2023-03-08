cwlVersion: v1.2
class: Operation
id: gunzip
label: gunzip

requirements:
  DockerRequirement:
    dockerPull: 'ubuntu:xenial'

inputs:
  input_fasta:
    type: File

outputs:
  fa:
    type: File
