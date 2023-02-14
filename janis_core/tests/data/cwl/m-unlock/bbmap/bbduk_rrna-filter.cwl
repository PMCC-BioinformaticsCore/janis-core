#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool

label: "Filter rRNA from reads"

doc: |
    Filter rRNA reads from paired end reads using BBmaps bbduk tool

requirements:
 - class: InlineJavascriptRequirement

hints:
  SoftwareRequirement:
    packages:
      bbmap:
        version: ["38.98"]
        specs: ["https://anaconda.org/bioconda/bbmap"]

  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/bbmap:38.98

stderr: rRNA-filter_summary.txt

inputs:
  identifier:
    type: string
    doc: Identifier for this dataset used in this workflow
    label: identifier used
  memory:
    type: int?
    default: 8
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
    type: File?
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
    default: "/riboKmers.fa.gz"

arguments:
  - prefix: "-Xmx"
    separate: false
    valueFrom: $(inputs.memory)M
  - prefix: "out="
    separate: false
    valueFrom: $(inputs.identifier)_rRNA-filt_1.fq.gz
  - prefix: "out2="
    separate: false
    valueFrom: $(inputs.identifier)_rRNA-filt_2.fq.gz
  - prefix: "stats="
    separate: false
    valueFrom: $(inputs.identifier)_rRNA-filter_stats.txt

baseCommand: [bbmap.sh]

outputs:
  out_forward_reads:
    type: File
    outputBinding:
      glob: $(inputs.identifier)_rRNA-filt_1.fq.gz
  out_reverse_reads:
    type: File
    outputBinding:
      glob: $(inputs.identifier)_rRNA-filt_2.fq.gz

  stats_file:
    type: File
    outputBinding:
      glob: $(inputs.identifier)_rRNA-filter_stats.txt

  summary:
    type: stderr


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
s:dateModified: "2022-02-22"
s:license: https://spdx.org/licenses/Apache-2.0 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"


$namespaces:
  s: https://schema.org/