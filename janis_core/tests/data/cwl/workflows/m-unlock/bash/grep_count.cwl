#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool
label: "Perform a grep count on a file"

hints:
  DockerRequirement:
    dockerPull: debian:buster
    
requirements:
  - class: InlineJavascriptRequirement
  - class: ShellCommandRequirement

baseCommand: []

stdout: "count.txt"

inputs:
  binary:
    type: string
    default: "grep"
    inputBinding:
      position: 2
      shellQuote: False
  count:
    type: string
    default: "-c"
    inputBinding:
      position: 3
      shellQuote: False
  grep:
    type: string
    inputBinding:
      position: 4
  infile:
    type: File
    inputBinding:
      shellQuote: False 
      position: 5
  pipe:
    type: string
    default: "||"
    inputBinding:
      shellQuote: False 
      position: 6
  status:
    type: string
    default: "true"
    inputBinding:
      shellQuote: False 
      position: 7

outputs:
  matches:
    type: int
    outputBinding:
      glob: count.txt
      loadContents: true
      outputEval: $(parseInt(self[0].contents))

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