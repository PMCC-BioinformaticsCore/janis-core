cwlVersion: v1.0
class: CommandLineTool
id: bedops_union
label: bedops_union

doc: |
  Concat BED files

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - $(inputs.input)

hints:
  DockerRequirement:
    dockerPull: 'quay.io/biocontainers/bedops:2.4.39--h7d875b9_1'

baseCommand: [bedops, -u]

inputs:
  input:
    type: 'File[]'
    inputBinding:
      position: 1

stdout: $(inputs.input[0].basename.split(".")[1]).all.filtered.bed

outputs:
  output:
    type: stdout
