cwlVersion: v1.2
class: Operation
id: codex
label: codex

doc: |
  CODEX2

requirements:
  DockerRequirement:
    dockerPull: 'migbro/codex2:3.8'

inputs:
  input:
    type: 'File[]'
  mapping:
    type: File
  bed:
    type: File
  chromosome:
    type: string
    
outputs:
  output:
    type: File