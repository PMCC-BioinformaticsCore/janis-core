#!/usr/bin/env cwltool

cwlVersion: v1.2
class: CommandLineTool

# baseCommand: ["bash", "script.sh"]
baseCommand: [carve]

label: "CarveMe"
doc: |
    CarveMe is a python-based tool for genome-scale metabolic model reconstruction.
    (Workflow will quit as successful even though no model can be created. Check messages.)
        
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entry: "$({class: 'Directory', listing: []})"
        entryname: "carve_output"
        writable: true
      - entryname: script.sh
        entry: |-
          #!/bin/bash
          source /root/miniconda/bin/activate
          conda init bash
          conda activate /unlock/infrastructure/conda/carveme/cplex/carveme_1.5.1
          carve $@

hints:
  SoftwareRequirement:
    packages:
      carveme:
        version: ["1.5.1"]
        specs: ["https://anaconda.org/bioconda/carveme"]


outputs:
  carveme_gem:
    label: CarveMe GEM
    doc: CarveMe metabolic model Output SBML in sbml-fbc2 format
    type: File?
    outputBinding:
      glob: $(inputs.protein_file.nameroot).GEM.xml

inputs:
  protein_file:
    type: File
    label: Input fasta file
    doc: Proteins sequence file in FASTA format. 
    inputBinding:
        position: 0

  gapfill:
    type: string?
    label: Gap fill
    doc: Gap fill model for given media
    inputBinding:
      prefix: --gapfill
    # default: "M8"

  mediadb:
    type: File?
    label: Media database
    doc: Media database file
    inputBinding:
      prefix: --mediadb

  init:
    type: File?
    label: Initial media
    doc: Initialize model with given medium
    inputBinding:
      prefix: --init


arguments:
  - "--fbc2"
  - prefix: "--output"
    valueFrom: $(inputs.protein_file.nameroot).GEM.xml

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
s:dateCreated: "2022-06-00"
s:license: https://spdx.org/licenses/Apache-2.0 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"

$namespaces:
  s: https://schema.org/
