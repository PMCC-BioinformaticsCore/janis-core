#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool

label: "Filter from reads"

doc: |
    Filter reads using BBmaps bbduk tool (paired-end only)

requirements:
 - class: InlineJavascriptRequirement

hints:
  SoftwareRequirement:
    packages:
      bbmap:
        version: ["39.01"]
        specs: ["https://anaconda.org/bioconda/bbmap"]

  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/bbmap:39.01

baseCommand: [ bbduk.sh ]

inputs:
  identifier:
    type: string
    doc: Identifier for this dataset used in this workflow
    label: identifier used
  memory:
    type: int?
    default: 8000
  threads:
    type: int?
    default: 1
    inputBinding:
      prefix: 'threads='
      separate: false
  forward_reads:
    type: File
    inputBinding:
      prefix: "in="
      separate: false
  reverse_reads:
    type: File
    inputBinding:
      prefix: "in2="
      separate: false
  kmersize:
    type: int
    inputBinding:
      prefix: "k="
      separate: false
    default: 31
  reference:
    doc: Reference contamination fasta file (can be compressed)
    label: Reference contamination file
    type: string?
    inputBinding:
      prefix: "ref="
      separate: false

stderr: $(inputs.identifier)_bbduk-summary.txt

arguments:
  - prefix: "-Xmx"
    separate: false
    valueFrom: $(inputs.memory)M
  - prefix: "out="
    separate: false
    valueFrom: $(inputs.identifier)_1.fq.gz
  - prefix: "out2="
    separate: false
    valueFrom: $(inputs.identifier)_2.fq.gz
  - prefix: "stats="
    separate: false
    valueFrom: $(inputs.identifier)_bbduk-stats.txt

outputs:
  out_forward_reads:
    type: File
    outputBinding:
      glob: $(inputs.identifier)_1.fq.gz
  out_reverse_reads:
    type: File
    outputBinding:
      glob: $(inputs.identifier)_2.fq.gz
  stats_file:
    type: File
    outputBinding:
      glob: $(inputs.identifier)_bbduk-stats.txt
  summary:
    type: File?
    outputBinding:
      glob: $(inputs.identifier)_bbduk-summary.txt 

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
s:dateModified: "2023-02-10"
s:license: https://spdx.org/licenses/Apache-2.0 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"


$namespaces:
  s: https://schema.org/