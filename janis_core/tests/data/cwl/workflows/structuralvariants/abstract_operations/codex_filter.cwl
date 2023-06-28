cwlVersion: v1.2
class: Operation
id: codex_filter
label: codex_filter

inputs:
  input:
    type: File
  samples:
    type: File
  min_len:
    type: string
  max_len:
    type: string
  min_lratio:
    type: string

outputs:
  output:
    type: File