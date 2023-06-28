cwlVersion: v1.0
class: CommandLineTool
id: bgzip
label: bgzip

doc: |
  Compressing BED file

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - $(inputs.input)

hints:
  DockerRequirement:
    dockerPull: 'quay.io/biocontainers/tabix:0.2.6--ha92aebf_0'

baseCommand: [bgzip, -c]

inputs:
  input:
    type: File
    inputBinding:
      position: 1

stdout: $(inputs.input.basename).gz

outputs:
  output:
    type: stdout
