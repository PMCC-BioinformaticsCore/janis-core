#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0
label: Extract Strelka Somatic AD DP
doc: |2-

   - Extract and calculate AD and AF value for each variant (both SNVs and INDELs)
  Based on https://github.com/Illumina/strelka/blob/v2.9.x/docs/userGuide/README.md#somatic
          

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: michaelfranklin/pmacutil:0.1.1

inputs:
- id: vcf
  label: vcf
  doc: input vcf
  type: File
  inputBinding:
    prefix: -i
- id: outputFilename
  label: outputFilename
  doc: output vcf
  type:
  - string
  - 'null'
  default: generated.vcf
  inputBinding:
    prefix: -o

outputs:
- id: out
  label: out
  type: File
  outputBinding:
    glob: generated.vcf
    loadContents: false
stdout: _stdout
stderr: _stderr

baseCommand: extract_strelka_somatic_DP_AF.py
arguments: []
id: extractStrelkaSomaticADDP
