#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.2

requirements:
- class: InlineJavascriptRequirement
- class: StepInputExpressionRequirement
- class: MultipleInputFeatureRequirement

inputs:
- id: inp1
  type: string

outputs: []

steps:
- id: stp1
  in:
  - id: inp
    source:
    - inp1
    linkMerge: merge_nested
  run: tools/ArrayStepTool.cwl
  out:
  - id: out
id: test_add_single_to_array_edge
