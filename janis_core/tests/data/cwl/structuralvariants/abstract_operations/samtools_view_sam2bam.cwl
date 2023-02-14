cwlVersion: v1.2
class: Operation
id: samtools_view_sam2bam
label: samtools_view_sam2bam

doc: |
  Convert SAM to BAM

requirements:
  DockerRequirement:
    dockerPull: 'quay.io/biocontainers/samtools:1.5--2'

inputs:
  input:
    type: File
  threads:
    type: 'int?'

outputs:
  output:
    type: File