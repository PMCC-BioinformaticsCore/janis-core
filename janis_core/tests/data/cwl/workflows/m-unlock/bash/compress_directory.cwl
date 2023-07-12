#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool
label: "Compress a directory (tar)"

hints:
  DockerRequirement:
    dockerPull: debian:buster
    
# tar -czf PROVENANCE.tar.gz -C $destination_path PROVENANCE
baseCommand: [tar, czfh]

arguments:
  - valueFrom: $(inputs.indir.basename).tar.gz
  - valueFrom: "-C"
  - valueFrom: $(inputs.indir.path)/..
  - valueFrom: $(inputs.indir.basename)

inputs:
  indir:
    type: Directory

outputs:
  outfile:
    type: File
    outputBinding:
        glob: $(inputs.indir.basename).tar.gz

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