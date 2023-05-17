#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

label: MaxBin2

doc: MaxBin2 is a software for binning assembled metagenomic sequences based on an Expectation-Maximization algorithm. 

requirements:
  - class: InlineJavascriptRequirement

hints:
  SoftwareRequirement:
    packages:
      maxbin2 :
        version: ["2.2.7"]
        specs: ["https://anaconda.org/bioconda/maxbin2"]
  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/maxbin2:2.2.7

baseCommand: [ "run_MaxBin.pl" ]

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
      prefix: -thread
  contigs:
    type: File
    doc: Input assembly in fasta format
    label: Input assembly
    inputBinding:
      prefix: -contig  
  abundances:
    type: File
    doc: Abundances file
    label: Abundances
    inputBinding:
      prefix: -abund

arguments:
  - prefix: "-out"
    valueFrom: $(inputs.identifier)_MaxBin2.bin

stdout: $(inputs.identifier)_MaxBin2.log

outputs:
  bins:
    type: File[]
    label: Bins
    doc: Bin fasta files. The XX bin. XX are numbers, e.g. out.001.fasta
    outputBinding:
      glob: "*.fasta"
  summary:
    type: File
    label: MaxBin2 Summary
    doc: Summary file describing which contigs are being classified into which bin.
    outputBinding:
      glob: $(inputs.identifier)_MaxBin2.bin.summary
  log:
    type: File
    label: Log
    doc: Log file recording the core steps of MaxBin algorithm
    outputBinding:
      glob: $(inputs.identifier)_MaxBin2.log
  markers:
    type: File
    label: Markers
    doc: Marker gene presence numbers for each bin. This table is ready to be plotted by R or other 3rd-party software.
    outputBinding:
      glob: $(inputs.identifier)_MaxBin2.bin.marker
  # noclass:
  #   type: File
  #   label: Not classfied seqs
  #   doc: All sequences that pass the minimum length threshold but are not classified successfully.
  #   outputBinding:
  #     glob: $(inputs.identifier)_MaxBin2.noclass
  # tooshort:
  #   type: File
  #   label: Too short
  #   doc: All sequences that do not meet the minimum length threshold.
  #   outputBinding:
  #     glob: $(inputs.identifier)_MaxBin2.tooshort

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
s:dateCreated: "2022-08-00"
s:license: https://spdx.org/licenses/Apache-2.0
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"

$namespaces:
  s: https://schema.org/
