cwlVersion: v1.2
class: Operation
id: bwa_mem_paired
label: bwa_mem_paired

requirements:
  DockerRequirement:
    dockerPull: 'quay.io/biocontainers/bwa:0.7.17--h84994c4_5'

inputs:
  reads:
    type: 'File[]'
  reference_genome:
    type: File
    secondaryFiles:
      - .fai
      - .amb
      - .ann
      - .bwt
      - .pac
      - .sa
  threads:
    type: 'int?'
  read_group:
    type: 'string?'

outputs:
  output:
    type: File