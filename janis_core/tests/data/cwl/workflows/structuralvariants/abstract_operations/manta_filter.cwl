cwlVersion: v1.2
class: Operation
id: manta_filter
label: manta_filter

inputs:
  input:
    type: File
  samples:
    type: File
  min_len:
    type: string
  max_len:
    type: string
  min_q:
    type: string

outputs:
  output:
    type: File