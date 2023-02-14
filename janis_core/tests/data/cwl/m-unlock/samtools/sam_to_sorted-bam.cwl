#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

label: "sam to sorted bam"
doc: |
  samtools view -@ $2 -hu $1 | samtools sort -@ $2 -o $3.bam

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entry: "$({class: 'Directory', listing: []})"
        entryname: "minimap_run"
        writable: true
      - entryname: script.sh
        entry: |-
          #!/bin/bash
          samtools view -@ $2 -hu $3 | samtools sort -@ $2 -o $1.sorted.bam

baseCommand: [ bash, -x, script.sh ]

hints:
  SoftwareRequirement:
    packages:
      samtools:
        version: ["1.15.1"]
        specs: ["https://anaconda.org/bioconda/samtools"]

  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/samtools:1.15.1

inputs:
  identifier:
    type: string
    doc: Identifier for this dataset used in this workflow
    label: identifier used
    inputBinding:
      position: 1

  sam:
    type: File
    doc: unsorted sam file
    label: unsorted sam file
    inputBinding:
      position: 3
  
  threads:
      type: int
      doc: number of cpu threads to use
      label: cpu threads
      default: 1
      inputBinding:
        position: 2


outputs:
  sortedbam: 
    type: File
    outputBinding:
      glob: $(inputs.identifier).sorted.bam       

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
