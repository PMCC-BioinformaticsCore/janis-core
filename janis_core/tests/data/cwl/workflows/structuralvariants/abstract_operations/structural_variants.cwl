cwlVersion: v1.2
class: Operation
id: structural_variants
label: structural_variants

doc: |
  Convert a VCF file to a BEDPE file

requirements:
  DockerRequirement:
    dockerPull: 'quay.io/biocontainers/bioconductor-structuralvariantannotation:1.6.0--r40_0'

inputs:
  input:
    type: File

outputs:
  output:
    type: File