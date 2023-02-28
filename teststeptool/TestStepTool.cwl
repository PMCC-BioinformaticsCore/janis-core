#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.2
label: TestStepTool

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: container/override

inputs:
- id: input1
  label: input1
  type: string
  inputBinding:
    position: 0
- id: input2
  label: input2
  type:
  - string
  - 'null'
  inputBinding:
    position: 1
- id: input3
  label: input3
  type:
  - string
  - 'null'
  inputBinding:
    position: 2
- id: input4
  label: input4
  type:
  - string
  - 'null'
  inputBinding:
    position: 3

outputs:
- id: out
  label: out
  type: string
  outputBinding:
    glob: '*'
    loadContents: false
stdout: _stdout
stderr: _stderr

baseCommand: echo
arguments: []

hints:
- class: ToolTimeLimit
  timelimit: |-
    $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
id: TestStepTool
