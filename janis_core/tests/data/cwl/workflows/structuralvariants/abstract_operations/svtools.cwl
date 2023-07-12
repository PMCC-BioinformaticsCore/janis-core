cwlVersion: v1.2
class: Operation
id: svtools
label: svtools

doc: |
  Convert a VCF file to a BEDPE file

requirements:
  DockerRequirement:
    dockerPull: 'quay.io/biocontainers/svtools:0.5.1--py_0'

inputs:
  input:
    type: File

outputs:
  output:
    type: File