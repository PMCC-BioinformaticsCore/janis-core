#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.2
label: EchoTestTool

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement

inputs:
- id: inp
  label: inp
  type: string
  inputBinding:
    position: 0

outputs:
- id: out
  label: out
  type: stdout
stdout: _stdout
stderr: _stderr

baseCommand: echo
arguments: []

hints:
- class: ToolTimeLimit
  timelimit: |-
    $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
id: EchoTestTool
