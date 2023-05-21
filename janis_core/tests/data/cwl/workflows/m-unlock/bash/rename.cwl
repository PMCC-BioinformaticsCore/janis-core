#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool
baseCommand: mv
label: Changes the name of a file and returns the renamed file as output

hints:
  DockerRequirement:
    dockerPull: debian:buster
    
requirements:
  InlineJavascriptRequirement: {}

inputs:
  newname:
    type: string
    inputBinding:
      position: 2
  rename_this:
    type: File
    inputBinding:
      position: 1

outputs:
  file_renamed:
    type: File
    outputBinding:
      glob: $(inputs.newname)