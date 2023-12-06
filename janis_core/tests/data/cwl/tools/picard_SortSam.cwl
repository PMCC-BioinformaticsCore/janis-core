#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/picard:2.22.2--0
  SoftwareRequirement:
    packages:
      picard:
        version: [ "2.22.2" ]
        specs: [ https://identifiers.org/biotools/picard_tools ]
requirements:
  InlineJavascriptRequirement: {}

baseCommand: [ picard, SortSam ]

arguments:
  - prefix: OUTPUT=
    separate: false
    valueFrom: |
      ${ if(inputs.sort_order == "coordinate") { return (inputs.alignments.nameroot)+".bam";} else { return (inputs.alignments.nameroot)+".sam"; } }

inputs:
  alignments:
    type: File
    inputBinding:
      prefix: INPUT=
      separate: false

  sort_order:
    type:
      - 'null'
      - type: enum
        symbols:
          - queryname
          - coordinate
          - duplicate
    default: coordinate
    doc: 'coordinate (bam) or queryname (sam)'
    inputBinding:
      prefix: SORT_ORDER=
      separate: false

  validation_stringency:
    default: LENIENT
    doc: Validation stringency for all SAM files read by this program.  Setting stringency
      to SILENT can improve performance when processing a BAM file in which variable-length
      data (read, qualities, tags) do not otherwise need to be decoded.
    type:
    - 'null'
    - type: enum
      symbols:
      - STRICT
      - LENIENT
      - SILENT
    inputBinding:
      prefix: VALIDATION_STRINGENCY=
      separate: false

outputs:
  sorted_alignments:
    type: File
    format: |-
      ${ if(inputs.sort_order == "coordinate") { return "http://edamontology.org/format_2572";} else { return "http://edamontology.org/format_2573"; } }
    outputBinding:
      glob: '*.*am'

$namespaces:
  edam: http://edamontology.org/
$schemas:
  - https://edamontology.org/EDAM_1.18.owl
