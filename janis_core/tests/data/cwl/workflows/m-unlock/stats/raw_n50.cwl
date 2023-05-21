#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool
label: Fasta statistics
doc: Fasta statistics like N50, total length, etc..

hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/idba:1.1.3--1

baseCommand: [raw_n50]

stdout: $(inputs.identifier)_stats.txt

inputs:
  identifier:
    type: string
    doc: Identifier for this dataset used in this workflow
    label: identifier used  
  input_fasta:
    type: File
    label: Input fasta
    doc: Input multi fasta file
    inputBinding:
        position: 1

outputs:
  output:
    type: File
    outputBinding:
        glob: $(inputs.identifier)_stats.txt

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
s:dateCreated: "2022-00-06"
s:license: https://spdx.org/licenses/Apache-2.0 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"

$namespaces:
  s: https://schema.org/
