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
        entryname: "fastq_to_fasta"
        writable: true
      - entryname: script.sh
        entry: |-
          #!/bin/bash
          zcat $@ | sed -n '1~4s/^@/>/p;2~4p'

baseCommand: [bash, script.sh]

stdout: $(inputs.identifier).fasta

inputs:
  identifier:
    type: string
    doc: Name of the output file
    label: output file name
  fastq_files:
      type: File[]
      doc: file list to be concatenated
      label: file list
      inputBinding:
        position: 0

outputs:
  fasta_out:
    type: stdout

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
