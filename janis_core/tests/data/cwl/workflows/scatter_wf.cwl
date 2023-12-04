

cwlVersion: v1.0
class: Workflow

requirements:
  ScatterFeatureRequirement: {}

inputs:
  message_array: string[]

steps:
  echo:
    run: ../tools/echo.cwl
    scatter: text
    in:
      text: message_array
    out: []

outputs: []
