#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

label: Filtlong

doc: |  
      Filtlong is a tool for filtering long reads by quality. It can take a set of long reads and produce a smaller, better subset. 
      It uses both read length (longer is better) and read identity (higher is better) when choosing which reads pass the filter.

hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/filtlong:0.2.1--hd03093a_1    

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entry: "$({class: 'Directory', listing: []})"
        entryname: "filtlong"
        writable: true
      - entryname: script.sh
        entry: |-
          #!/bin/bash
          outname=$1
          longreads=$2
          shift;shift;
          filtlong $longreads $@ 2> >(tee -a $outname.filtlong.log>&2) | gzip > $outname.fastq.gz

baseCommand: [ bash, -x, script.sh ]

outputs:
  output_reads:
    type: File?
    outputBinding:
      glob: $(inputs.output_filename).fastq.gz
  log:
    type: File?
    outputBinding:
      glob: $(inputs.output_filename).filtlong.log

inputs:
  output_filename:
    type: string
    label: Output filename
    doc: Output filename (fastq.gz will be added by default)
    inputBinding:
      position: 1

  long_reads:
    type: File
    label: Long reads
    doc: Long reads in fastq format
    inputBinding:
      position: 2

  target_bases:
    type: int?
    label: Target bases
    doc: Keep only the best reads up to this many total bases
    inputBinding:
      prefix: --target_bases
      position: 3
  keep_percent:
    type: float?
    label: Keep percentage
    doc: Keep only this percentage of the best reads (measured by bases)
    inputBinding:
      prefix: --keep_percent
      position: 4
  minimum_length:
    type: int?
    label: Minimum length
    doc: Minimum read length threshold
    inputBinding:
      prefix: --min_length
      position: 5
  maximum_length:
    type: int?
    label: Maximum length
    doc: Maximum read length threshold
    inputBinding:
      prefix: --max_length
      position: 6
  
  min_mean_q:
    type: float?
    label: Minimum mean quality
    doc: Minimum mean quality threshold
    inputBinding:
      prefix: --min_mean_q
      position: 7
  min_window_q:
    type: float?
    label: Minimum window quality
    doc: Minimum window quality threshold
    inputBinding:
      prefix: --min_window_q
      position: 8

  trim:
    type: boolean?
    label: Trim 
    doc: Trim non-k-mer-matching bases from start/end of reads
    inputBinding:
      prefix: --trim
      position: 9
  split:
    type: boolean?
    label: Split 
    doc: Split reads at this many (or more) consecutive non-k-mer-matching bases
    inputBinding:
      prefix: --trim
      position: 10

#  External references (if provided, read quality will be determined using these instead of from the Phred scores)
  forward_reads:
    type: File?
    label: Forward reads
    doc: Forward reference Illumina reads in FASTQ format
    inputBinding:
      prefix: -illumina_1
      position: 11
  reverse_reads:
    type: File?
    label: Reverse reads
    doc: Reverse reference Illumina reads in FASTQ format
    inputBinding:
      prefix: -illumina_2
      position: 12
  assembly:
    type: File?
    label: Reference assembly
    doc: Reference assembly in FASTA format
    inputBinding:
      prefix: --assembly
      position: 13

#  score weights (control the relative contribution of each score to the final read score):
  length_weight:
    type: float?
    label: Length weight
    doc: Weight given to the length score (default; 1)
    inputBinding:
      prefix: --length_weight
      position: 14
  mean_q_weight:
    type: float?
    label: Mean quality weight
    doc: Weight given to the mean quality score (default; 1)
    inputBinding:
      prefix: --mean_q_weight
      position: 15
  window_q_weight:
    type: float?
    label: Mean window weight
    doc: Weight given to the window quality score (default; 1)
    inputBinding:
      prefix: --window_q_weight
      position: 16

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
s:dateCreated: "2023-01-03"
s:license: https://spdx.org/licenses/Apache-2.0
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"

$namespaces:
  s: https://schema.org/
