#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.2
label: id

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement

inputs:
- id: inp
  label: inp
  type:
  - string
  - 'null'

outputs:
- id: out
  label: out
  type: string
  outputBinding:
    glob: $(inputs.inp)
    loadContents: false
stdout: _stdout
stderr: _stderr
arguments: []

hints:
- class: ToolTimeLimit
  timelimit: |-
    $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
id: id
