#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

label: "iRODS Sequence information"

hints:
  SoftwareRequirement:
    packages:
      python:
        version: ["3.10.6"]
        specs: ["https://anaconda.org/conda-forge/python"]
  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/python:3.10.6

requirements:
 - class: InlineJavascriptRequirement
 - class: InitialWorkDirRequirement
   listing:
      - entryname: ignore.txt
        entry: |-
          empty file

doc: |
    Processes sequence file for read information

inputs:
  fastq:
    type: string
    doc: Fastq file path in iRODS
    label: iRODs fastq file path
    inputBinding:
      prefix: -file
      position: 100
      

arguments: ["python3", "/unlock/infrastructure/scripts/reads.py"]

outputs:
  files:
    type: File
    outputBinding:
      glob: "ignore.txt"

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
s:license: https://spdx.org/licenses/Apache-2.0 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"


$namespaces:
  s: https://schema.org/
