#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

label: DAS Tool

doc: Recovery of genomes from metagenomes via a dereplication, aggregation and scoring strategy. 

requirements:
  - class: InlineJavascriptRequirement

hints:
  SoftwareRequirement:
    packages:
      dastool :
        version: ["1.1.5"]
        specs: ["https://anaconda.org/bioconda/das_tool"]
  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/das_tool:1.1.5

baseCommand: [ "DAS_Tool" ]

inputs:
  identifier:
    type: string
    doc: Identifier for this dataset used in this workflow
    label: identifier used
  threads:
    type: int?
    label: Number of threads to use
    default: 1
    inputBinding:
      position: 0
      prefix: --threads
  assembly:
    type: File
    doc: Input assembly in fasta format
    label: Input assembly
    inputBinding:
      prefix: --contigs
  bin_tables:
    type: File[]
    doc: Comma separated list of tab separated contigs to bin tables.
    label: Bin-Contig tables
    inputBinding:
      itemSeparator: ","
      prefix: --bins
  binner_labels:
    type: string
    doc: Comma separated list of the binning prediction tool names.
    label: Binner labels
    inputBinding:
      prefix: --labels

  write_bins:
    type: boolean?
    doc: Export bins as fasta files.
    label: Write bins
    inputBinding:
      prefix: --write_bins
    default: true
  write_unbinned:
    type: boolean?
    doc: Export unbinned contigs as fasta file
    label: Write unbinned
    inputBinding:
      prefix: --write_unbinned
    default: true

arguments:
  - prefix: "-o"
    valueFrom: $(inputs.identifier)

outputs:
  bin_dir:
    type: Directory?
    label: Bins
    doc: Bin fasta files.
    outputBinding:
      glob: $(inputs.identifier)_DASTool_bins
  summary:
    type: File
    label: DAS Tool run summary 
    doc: Summary
    outputBinding:
      glob: "*_DASTool_summary.tsv"
  contig2bin:
    type: File
    label:  Contig to bin
    doc: Contigs to bin file table
    outputBinding:
      glob: "*_DASTool_contig2bin.tsv"
  log:
    type: File
    label: Log
    doc: DASTool log file
    outputBinding:
      glob: "*_DASTool.log"

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
s:dateCreated: "2022-09-00"
s:license: https://spdx.org/licenses/Apache-2.0
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"

$namespaces:
  s: https://schema.org/
