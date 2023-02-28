#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.2
label: 'TEST: WorkflowWithStepInputExpression'

requirements:
- class: InlineJavascriptRequirement
- class: StepInputExpressionRequirement

inputs:
- id: mystring
  type:
  - string
  - 'null'
- id: mystring_backup
  type:
  - string
  - 'null'

outputs:
- id: out
  type: File
  outputSource: print/out

steps:
- id: print
  in:
  - id: _print_inp_mystring
    source: mystring
  - id: _print_inp_mystringbackup
    source: mystring_backup
  - id: inp
    valueFrom: |-
      $((inputs._print_inp_mystring != null) ? inputs._print_inp_mystring : inputs._print_inp_mystringbackup)
  run: tools/EchoTestTool_TEST.cwl
  out:
  - id: out
id: StepInputExpressionTestWF
