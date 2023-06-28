#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

doc: |
  NGtax2 output conversion to biom file

requirements:
  InlineJavascriptRequirement: {}

hints:
  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/scripts:1.0.1

inputs:
  job:
    type: File
    doc: Job file containing the folders to be processed
    label: Job file
  identifier:
    type: string
    doc: Prefix label for the output files
    label: Prefix label

 
baseCommand: ["python3", "/scripts/ngtax_to_biom.py"]

arguments:
  - "-j"
  - $(inputs.job)
  - "-i"
  - $(inputs.identifier)

outputs:
  biom_file:
    type: File
    outputBinding:
      glob: "*.biom"
  

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0001-8172-8981
    s:email: mailto:jasper.koehorst@wur.nl
    s:name: Jasper Koehorst

s:citation: https://m-unlock.nl
s:codeRepository: https://gitlab.com/m-unlock/cwl
s:dateCreated: "2021-00-00"
s:license: https://spdx.org/licenses/Apache-2.0 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"

$namespaces:
  s: https://schema.org/