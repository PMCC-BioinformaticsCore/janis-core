cwlVersion: v1.0
class: CommandLineTool
id: tabix
label: tabix

doc: |
  Indexing BED file

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - $(inputs.input)

hints:
  DockerRequirement:
    dockerPull: 'quay.io/biocontainers/tabix:0.2.6--ha92aebf_0'

baseCommand: [tabix]

inputs:
  input:
    type: File
    inputBinding:
      position: 1

outputs:
  output:
    type: File
    secondaryFiles: .tbi
    outputBinding:
      glob: $(inputs.input.basename)
