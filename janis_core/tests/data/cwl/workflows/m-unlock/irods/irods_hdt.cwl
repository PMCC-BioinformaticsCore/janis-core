#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

label: "iRODS HDT creator"

requirements:
 - class: InlineJavascriptRequirement
 - class: InitialWorkDirRequirement
   listing:
      - entryname: ignore.txt
        entry: |-
          empty file

doc: |
    Collects RDF (Turtle) files from given folder (recursively) and creates a HDT file

inputs:
  folder:
    type: string
    doc: Folder path in iRODS
    label: iRODs folder path
    inputBinding:
      position: 100

arguments: ["java","-jar", "/unlock/infrastructure/binaries/irods/RDFManager.jar","-hdt", "-folder"]

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
