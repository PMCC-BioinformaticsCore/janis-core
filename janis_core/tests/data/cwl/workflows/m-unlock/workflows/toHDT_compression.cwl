#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: Workflow
requirements:
   - class: StepInputExpressionRequirement
   - class: InlineJavascriptRequirement
   - class: MultipleInputFeatureRequirement

label: Conversion and compression of RDF files
doc: Workflow to convert a RDF file to the HDT format and GZIP compress it for long term storage

inputs:
  input:
    type: File
    doc: Genome sequence in EMBL format
    label: RDF input file
  output:
    type: string
    doc: Output filename
    label: HDT filename as output
  
  destination:
    type: string?
    label: Output Destination
    doc: Optional Output destination used for cwl-prov reporting.

steps:
############################
  tohdt:
    run: ../sapp/toHDT.cwl
    in:
      output: output
      input: input
    out: [output]
  #########################################
  #
  compress:
    run: ../bash/compress.cwl
    in:
      infile: tohdt/output
    out: [output]

outputs:
  output:
    type: File
    outputSource: compress/output
  
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
  