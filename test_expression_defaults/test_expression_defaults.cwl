#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.2

requirements:
- class: InlineJavascriptRequirement
- class: StepInputExpressionRequirement

inputs:
- id: inp
  type:
  - string
  - 'null'

outputs:
- id: out
  type: File
  outputSource: echo/out

steps:
- id: echo
  in:
  - id: _echo_inp_inp
    source: inp
  - id: inp
    valueFrom: |-
      $(("Hello, " + (inputs._echo_inp_inp != null) ? inputs._echo_inp_inp : ", Michael!"))
  run: tools/EchoTestTool_TEST.cwl
  out:
  - id: out
id: test_expression_defaults
