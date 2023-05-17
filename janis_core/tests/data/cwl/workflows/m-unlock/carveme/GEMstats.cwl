#!/usr/bin/env cwltool

cwlVersion: v1.2
class: CommandLineTool

baseCommand: ["bash","-x", "script.sh"]

label: "CarveMe GEMstats"
doc: |
    Small summary of a list of CarveMe genome-scale metabolic models in sbml-fbc2 format
    Contains; number of mets,reactions and genes

requirements:
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
      - entry: "$({class: 'Directory', listing: []})"
        entryname: "gemstats_output"
        writable: true
      - entryname: script.sh
        entry: |-
          #!/bin/bash
          identifier=$1
          shift;
          echo "Model Mets Reactions Genes" > $identifier\_CarveMe_GEMstats.tsv
          for file in "$@"
          do
            bash /unlock/infrastructure/scripts/GEMstats.sh $file
          done >> $identifier\_CarveMe_GEMstats.tsv

inputs:
  identifier:
    type: string
    doc: Identifier for this dataset used in this workflow
    label: identifier used
    inputBinding:
      position: 1
  carveme_gems:
    type: File[]
    label: CarveMe GEMs
    doc: List of CarveMe metabolic models in sbml-fbc2 format.
    inputBinding:
      position: 2

outputs:
  carveme_GEMstats: 
    type: File
    outputBinding:
      glob: $(inputs.identifier)_CarveMe_GEMstats.tsv


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
s:dateCreated: "2022-01-01"
s:license: https://spdx.org/licenses/Apache-2.0 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"

$namespaces:
  s: https://schema.org/