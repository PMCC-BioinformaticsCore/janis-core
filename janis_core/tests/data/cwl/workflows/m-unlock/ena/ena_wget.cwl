#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool

label: "ENA retrieval"

requirements:
 - class: InlineJavascriptRequirement

baseCommand: []

arguments:
  - valueFrom: wget
  - valueFrom: https://www.ebi.ac.uk/ena/browser/api/embl/$(inputs.identifier)?download=true&gzip=true

inputs:
  identifier:
    type: string
    doc: Genome Accession number according to ENA
    label: Genome accession number (GCA_....)

outputs:
  gz:
    type: File?
    outputBinding:
      glob: $(inputs.identifier)*
      outputEval: ${self[0].basename=inputs.identifier + ".embl.gz"; return self;}


$namespaces:
 s: http://schema.org/