cwlVersion: v1.2
class: Operation
id: fastqc
label: fastqc

doc: |
  Run fastqc on raw reads in FASTQ format (single or paired end)

requirements:
  DockerRequirement:
    dockerPull: 'quay.io/biocontainers/fastqc:0.11.8--1'

inputs:
  fastq1:
    type: 'File[]'
  fastq2:
    type: 'File[]'
  threads:
    type: 'int?'

outputs:
  fastqc_zip:
    type:
      type: array
      items: File
  fastqc_html:
    type:
      type: array
      items: File