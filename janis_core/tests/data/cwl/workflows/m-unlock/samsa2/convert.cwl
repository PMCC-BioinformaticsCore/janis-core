#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

doc: |
  Samsa2 conversion workflow

requirements:
 - class: InlineJavascriptRequirement

inputs:
  inputfile:
    type: File
    doc: diamond refseq or seed result table with salltitles
    label: diamond tabular file
  identifier:
    type: string
    doc: Identifier for this dataset used in this workflow
    label: identifier used

baseCommand: [python3]

arguments:
  - "/unlock/infrastructure/scripts/samsa2/manager.py"
  - $(inputs.inputfile.path)
  - $(inputs.identifier)

outputs:
  output_seed:
    type: File?
    outputBinding:
      glob: "*.hierarchy"
  output_seed_reduced:
    type: File?
    outputBinding:
      glob: "*.reduced"
  output_refseq_function:
    type: File?
    outputBinding:
      glob: "*_function.tsv"
  output_refseq_organism:
    type: File?
    outputBinding:
      glob: "*_organism.tsv"

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
