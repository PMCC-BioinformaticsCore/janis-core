cwlVersion: v1.0
class: CommandLineTool
id: picard_markduplicates
label: picard_markduplicates

doc: |
  Removal of duplicates from aligned reads

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - $(inputs.input)

hints:
  DockerRequirement:
    dockerPull: 'quay.io/biocontainers/picard:2.10.6--py35_0'

baseCommand: [picard, MarkDuplicates]

arguments:
  - position: 2
    valueFrom: OUTPUT=$(inputs.input.nameroot).dedup.bam
  - position: 3
    valueFrom: METRICS_FILE=$(inputs.input.nameroot).dedup.metrics.txt
  - position: 4
    valueFrom: ASSUME_SORTED=TRUE

inputs:
  input:
    type: File
    inputBinding:
      position: 1
      prefix: 'INPUT='
      separate: false

outputs:
  alignments:
    type: File
    outputBinding:
      glob: '*.sorted.dedup.bam'
  metrics:
    type: File
    outputBinding:
      glob: '*.sorted.dedup.metrics.txt'
