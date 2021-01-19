#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.2
label: Hello, World!
doc: |
  This is the 'Hello, world' equivalent workflow that uses the Echo unix
  tool to log "Hello, World!" to the console, and collects the result.

  This is designed to be the first example that you can run with janis, ie:
      
  ``janis run hello``

requirements:
- class: InlineJavascriptRequirement
- class: StepInputExpressionRequirement

inputs:
- id: inp
  type: string
  default: Hello, world!

outputs:
- id: out
  type: File
  outputSource: hello/out

steps:
- id: hello
  label: Echo
  in:
  - id: inp
    source: inp
  run: tools/echo_v1_0_0.cwl
  out: ["out"]
id: hello
