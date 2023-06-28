
cwlVersion: v1.0
class: CommandLineTool
id: svtools
label: svtools

doc: |
  Convert a VCF file to a BEDPE file

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - $(inputs.input)

hints:
  DockerRequirement:
    dockerPull: 'quay.io/biocontainers/svtools:0.5.1--py_0'

baseCommand: [svtools, vcftobedpe]

arguments:
  - position: 2
    prefix: '-o'
    valueFrom: $(inputs.input.nameroot.replace('vcf', 'bed'))

inputs:
  input:
    type: File
    inputBinding:
      position: 1
      prefix: '-i'

outputs:
  output:
    type: File
    outputBinding:
      glob: '*.bed'
