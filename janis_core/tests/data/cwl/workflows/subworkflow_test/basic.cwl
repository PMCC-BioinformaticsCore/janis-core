#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.2
label: basic tool for testing

requirements:
    - class: ShellCommandRequirement
    - class: InlineJavascriptRequirement
    - class: DockerRequirement
      dockerPull: ubuntu:latest

baseCommand: echo

arguments:
    - position: 0
      valueFrom: test:\\t:escaped:\\n:characters

inputs:
    inFile:
        type: File
        inputBinding:
            position: 1
    inString:
        type: string
        inputBinding:
            position: 2
    inSecondary:
        type: File
        secondaryFiles: [.fai, ^.dict, .amb, .ann, .bwt, .pac, .sa]
        inputBinding:
            position: 3
    inFileOptional:
        type: File?
        inputBinding:
            position: 4
    inStringOptional:
        type: string?
        inputBinding:
            position: 5
    inIntOptional:
        type: int?
        inputBinding:
            position: 6
    inIntOptional2:
        type: int?
        default: 10
        inputBinding:
            position: 6
outputs:
    out_stdout:
        type: stdout     




