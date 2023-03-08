cwlVersion: v1.2
class: Operation
id: manta
label: manta

doc: |
  Manta Structural Variant Caller

requirements:
  DockerRequirement:
    dockerPull: 'quay.io/biocontainers/manta:1.6.0--py27_0'

inputs:
  input:
    type: File
    secondaryFiles:
      - .bai
  reference_genome:
    type: File
    secondaryFiles:
      - .fai
      - .amb
      - .ann
      - .bwt
      - .pac
      - .sa
  exome:
    type: boolean
  regions:
    type: File
    secondaryFiles:
      - .tbi
  runDir:
    type: 'string?'

outputs:
  python:
    type: File
  output:
    type: File