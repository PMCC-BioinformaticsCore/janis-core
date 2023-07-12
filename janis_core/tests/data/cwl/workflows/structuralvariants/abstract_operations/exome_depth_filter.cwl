cwlVersion: v1.2
class: Operation
id: exomedepth_filter
label: exomedepth_filter

inputs:
  input:
    type: File
  samples:
    type: File
  min_len:
    type: string
  max_len:
    type: string
  min_bf:
    type: string

outputs:
  output:
    type: File