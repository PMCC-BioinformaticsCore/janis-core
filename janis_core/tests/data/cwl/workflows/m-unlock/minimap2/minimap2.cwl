#!/usr/bin/env cwltool

cwlVersion: v1.2
class: CommandLineTool

baseCommand: [minimap2]

label: "Minimap2"
doc: |
  "Minimap2 is a versatile sequence alignment program that aligns DNA or mRNA sequences against a large reference database. 
   Developed for mapping long reads from PacBio or Oxford Nanopore"

requirements:
  InlineJavascriptRequirement: {}

hints:
  SoftwareRequirement:
    packages:
      minimap2:
        version: ["2.24"]
        specs: ["https://anaconda.org/bioconda/minimap2"]
  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/minimap2:2.24

arguments:
  - valueFrom: "-a"
    position: 1
  - valueFrom: $(inputs.identifier).sam
    prefix: -o
    position: 2

inputs:
  threads:
    label: Number of CPU-threads used by minimap2.
    type: int?
    inputBinding:
      position: 3
      prefix: -t
  identifier:
    type: string
    doc: Identifier for this dataset used in this workflow
    label: identifier used
  reference:
    label: Target sequence in FASTQ/FASTA format (can be gzipped).
    type: File
    inputBinding:
      position: 4
  reads:
    label: Query sequence in FASTQ/FASTA format (can be gzipped).
    type: File[]
    inputBinding:
      position: 5
  preset:
    label: read type
    doc: |
      - map-pb/map-ont - PacBio CLR/Nanopore vs reference mapping
      - map-hifi - PacBio HiFi reads vs reference mapping
      - ava-pb/ava-ont - PacBio/Nanopore read overlap
      - asm5/asm10/asm20 - asm-to-ref mapping, for ~0.1/1/5% sequence divergence
      - splice/splice:hq - long-read/Pacbio-CCS spliced alignment
      - sr - genomic short-read mapping
    type: string
    inputBinding:
      position: 0
      prefix: -x

stderr: $(inputs.identifier)_minimap2.log

outputs:
  sam:
    label: Alignment of target vs query in SAM format
    type: File
    outputBinding:
      glob: $(inputs.identifier).sam
  log:
    label: minimap2 log
    type: File
    outputBinding:
      glob: $(inputs.identifier)_minimap2.log

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-5516-8391
    s:email: mailto:german.royvalgarcia@wur.nl
    s:name: Germán Royval
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
s:dateCreated: "2021-11-25"
s:license: https://spdx.org/licenses/Apache-2.0 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"

$namespaces:
  s: https://schema.org/
