#!/usr/bin/env cwl-runner

# TO BE REMOVED IN FAVOR OF PIGZ.CWL

cwlVersion: v1.2
class: CommandLineTool
label: "compress any file with gzip"

hints:
  DockerRequirement:
    dockerPull: debian:buster


baseCommand: [gzip, -c]

arguments:
  - valueFrom: $(inputs.infile)

stdout: $(inputs.infile.basename).gz

inputs:
  infile:
    type: File

outputs:
  output:
    type: File
    outputBinding:
        glob: $(inputs.infile.basename).gz

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