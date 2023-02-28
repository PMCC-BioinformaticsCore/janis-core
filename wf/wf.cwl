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
  type: string
  outputSource: _evaluate-output-out/out

steps:
- id: _evaluate-output-out
  in:
  - id: _inp
    source: inp
  run:
    class: ExpressionTool

    inputs:
    - id: _inp
      type:
      - string
      - 'null'
      loadContents: false

    outputs:
    - id: out
      type: string
    expression: '${return {out: inputs._inp }}'
  out:
  - out
id: wf
