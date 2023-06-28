#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool
label: "compress a file multithreaded with pigz"

hints:
  SoftwareRequirement:
    packages:
      pigz:
        version: ["2.6"]
        specs: ["https://anaconda.org/conda-forge/pigz"]
  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/pigz:2.6

baseCommand: [pigz, -c]

arguments:
  - valueFrom: $(inputs.inputfile)

stdout: $(inputs.inputfile.basename).gz

inputs:
  inputfile:
    type: File
  threads:
    type: int?
    default: 1
    inputBinding:
      prefix: '-p'

outputs:
  outfile:
    type: File
    outputBinding:
        glob: $(inputs.inputfile.basename).gz