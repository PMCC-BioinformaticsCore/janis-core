#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: Workflow
requirements:
   - class: StepInputExpressionRequirement
   - class: InlineJavascriptRequirement
   - class: MultipleInputFeatureRequirement

inputs:
  forward_primer:
    type: string
    doc: Forward primer used
    label: The forward primer used
  reverse_primer:
    type: string
    doc: Reverse primer used
    label: The reverse primer used
  forward_reads:
    type: File
    doc: Forward library file
    label: The forward library used
  reverse_reads: 
    type: File
    doc: Reverse library file
    label: The reverse library used
  mapping_file:
    type: File
    doc: Mapping file containing barcode information
    label: The mapping file
  destination:
    type: string?
    label: Output Destination
    doc: Optional Output destination used for cwl-prov reporting.
    
steps:
############################
  ngtax:
    run: ../ngtax/ngtax_demultiplexing.cwl
    in:
      forward_primer: forward_primer
      reverse_primer: reverse_primer
      forward_reads: forward_reads
      reverse_reads: reverse_reads
      mapping_file: mapping_file
    out: [output]


outputs:
  demultiplex_output:
    type: File[]
    outputSource: ngtax/output


s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0001-8172-8981
    s:email: mailto:jasper.koehorst@wur.nl
    s:name: Jasper Koehorst
  - class: s:Person
    s:identifier: https://orcid.org/0000-0001-9524-5964
    s:email: mailto:bart.nijsse@wur.nl
    s:name: Bart Nijsse

s:citation: https://m-unlock.nl
s:codeRepository: https://gitlab.com/m-unlock/cwl
s:dateCreated: "2020-00-00"
s:dateModified: "2022-05-00"
s:license: https://spdx.org/licenses/Apache-2.0 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"


$namespaces:
  s: https://schema.org/