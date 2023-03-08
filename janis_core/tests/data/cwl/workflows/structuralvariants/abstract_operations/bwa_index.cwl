cwlVersion: v1.2
class: Operation
id: bwa_index
label: bwa_index

requirements:
  DockerRequirement:
    dockerPull: 'quay.io/biocontainers/bwa:0.7.17--h84994c4_5'

inputs:
  generate_bwa_indexes:
    type: 'boolean?'
  input_fasta:
    type: File
  input_amb:
    type: 'File?'
  input_ann:
    type: 'File?'
  input_bwt:
    type: 'File?'
  input_pac:
    type: 'File?'
  input_sa:
    type: 'File?'
  algoType:
    type: 'string?'

outputs:
  amb:
    type: 'File?'
  ann:
    type: 'File?'
  bwt:
    type: 'File?'
  pac:
    type: 'File?'
  sa:
    type: 'File?'
