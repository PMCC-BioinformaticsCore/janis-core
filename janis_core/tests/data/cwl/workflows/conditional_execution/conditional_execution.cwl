#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: Workflow

label: |
  map medium and long reads (> 100 bp) against reference genome

inputs:
  in_reads:
    type: File
  in_reference: 
    type: File
  in_metadata:
    type: File?
  
outputs:
  out_alignments:
    type: File
    format: edam:format_2572  # BAM
    outputSource: align/aligned_reads
    label: "Alignments of the reads to the references genome"

  out_metadata:
    type: File
    outputSource: [ extract_metadata/output_metadata, in_metadata ]
    pickValue: the_only_non_null

steps:
  extract_metadata:
    run: extract_metadata.cwl
    label: "extract metadata from reads / reference"
    in:
      reads: in_reads
      reference: in_reference
      metadata: in_metadata
    when: $(inputs.metadata === null)
    out: [ output_metadata ]
  
  align: 
    run: align.cwl 
    in:
      reads: in_reads
      reference: in_reference
      metadata:
        source: [ extract_metadata/output_metadata, in_metadata ]
        pickValue: first_non_null
    out: [ aligned_reads ]

requirements:
  MultipleInputFeatureRequirement: {}
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}
#   SchemaDefRequirement:
#     types:
#       - $import: ReadGroupType.yml
 
$namespaces:
  edam: http://edamontology.org/
$schemas:
  - https://edamontology.org/EDAM_1.18.owl
