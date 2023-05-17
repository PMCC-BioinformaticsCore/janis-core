#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

label: "fastq to fasta"

doc: |
    Convert fastq file(s) to fasta format.
    zcat $@ | sed -n '1~4s/^@/>/p;2~4p'


requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entry: "$({class: 'Directory', listing: []})"
        entryname: "prepare_fasta_db"
        writable: true
      - entryname: script.sh
        entry: |-
          #!/bin/bash
          #zcat $@ | sed -n '1~4s/^@/>/p;2~4p'
          echo "bla"
          # echo $2
          # for f in $1/bin.*; do mv $f $2"_SemiBin_$f; done"

baseCommand: [bash, script.sh]

stdout: $(inputs.identifier).fasta

inputs:
  bin_directory:
      type: Directory
      doc: Directory with bins
      label: file list
      inputBinding:
        position: 1
  identifier:
    type: string
    doc: Name of the output file
    label: output file name
    inputBinding:
      position: 2

outputs:
  fasta_out:
    type: Directory?


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


