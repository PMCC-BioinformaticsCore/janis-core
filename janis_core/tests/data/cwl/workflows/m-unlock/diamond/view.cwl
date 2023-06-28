#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

doc: |
  Diamond workflow implementation

requirements:
 - class: InlineJavascriptRequirement

hints:
  SoftwareRequirement:
    packages:
      diamond :
        version: ["2.0.15"]
        specs: ["https://anaconda.org/bioconda/diamond"]
  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/diamond:2.0.15

inputs:
  inputfile:
    type: File
    doc: Diamond binary result file
    label: input file
    inputBinding:
        prefix: --daa

baseCommand: [ diamond ]

arguments:
  - "view"
  - "--outfmt"
  - "6"
  - "qseqid"
  - "sseqid"
  - "pident"
  - "length"
  - "mismatch"
  - "gapopen"
  - "qstart"
  - "qend"
  - "sstart"
  - "send"
  - "evalue"
  - "bitscore"
  - "stitle"
  - "--out"
  - "$(inputs.inputfile.basename).tsv"


outputs:
  output_diamond_tabular:
    type: File
    outputBinding:
      glob: $(inputs.inputfile.basename).tsv

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
