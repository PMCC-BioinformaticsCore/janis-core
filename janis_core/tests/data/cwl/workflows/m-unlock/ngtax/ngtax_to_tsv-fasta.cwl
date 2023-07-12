#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

doc: |
  NGtax2 output conversion to prepare for biom file and ASV fasta file

requirements:
  InlineJavascriptRequirement: {}

hints:
  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/scripts:1.0.1

inputs:
  input:
    type: File
    doc: NGTax2 turtle output file
    label: NGTax2 turtle file
  metadata:
    type: File?
    doc: UNLOCK assay metadata file
    label: Metadata file
  identifier:
    type: string
    doc: Identifier for this dataset used in this workflow
    label: identifier used
  fragment:
    type: string
    doc: Fragment identifier that was used
    label: fragment identifier

baseCommand: ["python3", "/scripts/ngtax_to_tsv-fasta.py"]

arguments:
  - "-t"
  - $(inputs.input.path)
  - "-i"
  - $(inputs.identifier)
  - "-f"
  - $(inputs.fragment)  
  - "-m"
  - $(inputs.metadata)  

outputs:
  picrust_tsv:
    type: File
    outputBinding:
      glob: "*.picrust.tsv"
  picrust_fasta:
    type: File
    outputBinding:
      glob: "*.picrust.fasta"
  physeq_asv:
    type: File
    outputBinding:
      glob: "*_asv.tsv"
  physeq_tax:
    type: File
    outputBinding:
      glob: "*_tax.tsv"
  physeq_seq:
    type: File
    outputBinding:
      glob: "*_seq.tsv"
  physeq_met:
    type: File
    outputBinding:
      glob: "*_met.tsv"

  

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
s:dateCreated: "2021-00-00"
s:license: https://spdx.org/licenses/Apache-2.0 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"

$namespaces:
  s: https://schema.org/