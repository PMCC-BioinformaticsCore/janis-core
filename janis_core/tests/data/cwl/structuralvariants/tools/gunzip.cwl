cwlVersion: v1.0
class: CommandLineTool
id: gunzip
label: gunzip

requirements:
  InlineJavascriptRequirement: {}
  ShellCommandRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - $(inputs.input_fasta)

hints:
  DockerRequirement:
    dockerPull: 'ubuntu:xenial'

baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      $(inputs.input_fasta.nameext == '.gz' ? 'gunzip -c' : 'echo gunzip')

inputs:
  input_fasta:
    type: File
    inputBinding:
      position: 2

stdout: $(inputs.input_fasta.basename.split('.')[0]).fa

outputs:
  fa:
    type: stdout
