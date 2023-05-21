

#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.2
label: basic tool for testing optional inputs

requirements:
    - class: DockerRequirement
      dockerPull: ubuntu:latest

baseCommand: echo

inputs:
    inFile:
        type: File?
        inputBinding:
            position: 1
    inFileArr:
        type: File[]?
        inputBinding:
            position: 2
    inSecondary:
        type: File?
        secondaryFiles: [.fai, ^.dict, .amb, .ann, .bwt, .pac, .sa]
        inputBinding:
            position: 3
    inString:
        type: string?
        inputBinding:
            position: 4
    inStringArr:
        type: string[]?
        inputBinding:
            position: 5
    inInt:
        type: int?
        inputBinding:
            position: 6
outputs:
    out_stdout:
        type: stdout     




