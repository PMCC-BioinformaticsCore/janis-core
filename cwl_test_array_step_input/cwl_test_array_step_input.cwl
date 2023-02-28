#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.2

requirements:
- class: InlineJavascriptRequirement
- class: StepInputExpressionRequirement
- class: MultipleInputFeatureRequirement

inputs:
- id: inp1
  type:
  - string
  - 'null'
- id: inp2
  type:
  - string
  - 'null'

outputs:
- id: out
  type:
    type: array
    items: File
  outputSource: print/out

steps:
- id: print
  in:
  - id: _print_inp_inp1
    source: inp1
  - id: _print_inp_inp2
    source: inp2
  - id: inp
    valueFrom: |-
      $([(inputs._print_inp_inp1 != null) ? inputs._print_inp_inp1 : "default1", (inputs._print_inp_inp2 != null) ? (inputs._print_inp_inp2 + "_suffix") : ""])
  run: tools/ArrayStepTool.cwl
  out:
  - id: out
id: cwl_test_array_step_input
