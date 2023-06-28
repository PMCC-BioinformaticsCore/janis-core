#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0
label: 'GATK4: LearnReadOrientationModel'
doc: TBD

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: broadinstitute/gatk:4.1.4.0

inputs:
- id: javaOptions
  label: javaOptions
  type:
  - type: array
    items: string
  - 'null'
- id: compression_level
  label: compression_level
  doc: |-
    Compression level for all compressed files created (e.g. BAM and VCF). Default value: 2.
  type:
  - int
  - 'null'
- id: f1r2CountsFiles
  label: f1r2CountsFiles
  doc: Counts for the read orientation of fragments
  type:
    type: array
    inputBinding:
      prefix: -I
    items: File
  inputBinding:
    position: 0
- id: numEmIterations
  label: numEmIterations
  doc: Amount of iterations for the em process before it bails
  type: int
  default: 30
  inputBinding:
    prefix: --num-em-iterations
    position: 1
- id: modelFileOut
  label: modelFileOut
  type:
  - string
  - 'null'
  default: generated.tar.gz
  inputBinding:
    prefix: -O
    position: 3

outputs:
- id: out
  label: out
  doc: Model file containing information about fragment orientations
  type: File
  outputBinding:
    glob: generated.tar.gz
    loadContents: false
stdout: _stdout
stderr: _stderr

baseCommand:
- gatk
- LearnReadOrientationModel
arguments:
- prefix: --java-options
  position: -1
  valueFrom: |-
    $("-Xmx{memory}G {compression} {otherargs}".replace(/\{memory\}/g, (([inputs.runtime_memory, 32, 4].filter(function (inner) { return inner != null })[0] * 3) / 4)).replace(/\{compression\}/g, (inputs.compression_level != null) ? ("-Dsamjdk.compress_level=" + inputs.compression_level) : "").replace(/\{otherargs\}/g, [inputs.javaOptions, []].filter(function (inner) { return inner != null })[0].join(" ")))
id: Gatk4LearnReadOrientationModel
