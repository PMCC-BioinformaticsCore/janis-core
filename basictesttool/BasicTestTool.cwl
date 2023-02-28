#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.2
label: Tool for testing translation

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
- class: EnvVarRequirement
  envDef:
  - envName: test1
    envValue: $(inputs.testtool)
- class: DockerRequirement
  dockerPull: ubuntu:latest

inputs:
- id: testtool
  label: testtool
  type: string
- id: arrayInp
  label: arrayInp
  type:
  - type: array
    items: string
  - 'null'

outputs:
- id: std
  label: std
  type: stdout
stdout: _stdout
stderr: _stderr

baseCommand: echo
arguments:
- position: 0
  valueFrom: test:\\t:escaped:\\n:characters\"

hints:
- class: ToolTimeLimit
  timelimit: |-
    $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
id: BasicTestTool
