#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.2
label: ArrayStepTool

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement

inputs:
- id: inp
  label: inp
  type:
    type: array
    items: string
  inputBinding:
    position: 1

outputs:
- id: out
  label: out
  type:
    type: array
    items: File
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
id: ArrayStepTool
