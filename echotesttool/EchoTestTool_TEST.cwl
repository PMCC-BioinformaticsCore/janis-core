#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.2
label: EchoTestTool

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: ubuntu:latest
- class: ResourceRequirement
  coresMin: |-
    $([inputs.runtime_cpu, (2 * inputs.outputFiles), 1].filter(function (inner) { return inner != null })[0])
  outdirMin: |-
    $([inputs.runtime_disk, 20].filter(function (inner) { return inner != null })[0])
  ramMin: |-
    $(Math.round((953.674 * [inputs.runtime_memory, ((inputs.inputFile.size / 1048576) > 1024) ? 4 : 2, 4].filter(function (inner) { return inner != null })[0])))

inputs:
- id: inputFile
  label: inputFile
  type: File
  inputBinding:
    position: 1
- id: outputFiles
  label: outputFiles
  type: int
- id: runtime_memory
  type: float?
- id: runtime_cpu
  type: int?
- id: runtime_disk
  type: int?
- id: runtime_seconds
  type: int?

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
    $([inputs.runtime_seconds, 60, 86400].filter(function (inner) { return inner != null })[0])
id: EchoTestTool
