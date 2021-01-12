#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.2
label: Echo
doc: |-
  The echo utility writes any specified operands, separated by single blank (` ') characters and followed by a newline (`
  ') character, to the standard output.

  Some shells may provide a builtin echo command which is similar or identical to this utility. Most notably, the builtin echo in sh(1) does not accept the -n option. Consult the builtin(1) manual page.

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: ubuntu@sha256:1d7b639619bdca2d008eca2d5293e3c43ff84cbee597ff76de3b7a7de3e84956

inputs:
- id: inp
  label: inp
  type: string
  inputBinding:
    position: 1
- id: include_newline
  label: include_newline
  doc: |-
    Do not print the trailing newline character.  This may also be achieved by appending `\c' to the end of the string, as is done by iBCS2 compatible systems.  Note that this option as well as the effect of `\c' are implementation-defined in IEEE Std 1003.1-2001 (``POSIX.1'') as amended by Cor. 1-2002.  Applications aiming for maximum portability are strongly encouraged to use printf(1) to suppress the newline character.
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: -n

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
id: echo
