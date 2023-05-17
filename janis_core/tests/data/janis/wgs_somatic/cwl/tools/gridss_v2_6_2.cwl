#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0
label: Gridss
doc: |
  GRIDSS: the Genomic Rearrangement IDentification Software Suite

  GRIDSS is a module software suite containing tools useful for the detection of genomic rearrangements.
  GRIDSS includes a genome-wide break-end assembler, as well as a structural variation caller for Illumina
  sequencing data. GRIDSS calls variants based on alignment-guided positional de Bruijn graph genome-wide
  break-end assembly, split read, and read pair evidence.

  GRIDSS makes extensive use of the standard tags defined by SAM specifications. Due to the modular design,
  any step (such as split read identification) can be replaced by another implementation that also outputs
  using the standard tags. It is hoped that GRIDSS can serve as an exemplar modular structural variant
  pipeline designed for interoperability with other tools.

  If you have any trouble running GRIDSS, please raise an issue using the Issues tab above. Based on feedback
  from users, a user guide will be produced outlining common workflows, pitfalls, and use cases.

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: gridss/gridss:2.6.2

inputs:
- id: bams
  label: bams
  type:
    type: array
    items: File
  inputBinding:
    position: 10
- id: reference
  label: reference
  type: File
  secondaryFiles:
  - .fai
  - .amb
  - .ann
  - .bwt
  - .pac
  - .sa
  - ^.dict
  inputBinding:
    prefix: --reference
    position: 1
- id: outputFilename
  label: outputFilename
  type:
  - string
  - 'null'
  default: generated.svs.vcf
  inputBinding:
    prefix: --output
    position: 2
- id: assemblyFilename
  label: assemblyFilename
  type:
  - string
  - 'null'
  default: generated.assembled.bam
  inputBinding:
    prefix: --assembly
    position: 3
- id: threads
  label: threads
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --threads
    valueFrom: |-
      $([inputs.runtime_cpu, 8, 1].filter(function (inner) { return inner != null })[0])
- id: blacklist
  label: blacklist
  type:
  - File
  - 'null'
  inputBinding:
    prefix: --blacklist
    position: 4
- id: tmpdir
  label: tmpdir
  type: string
  default: ./TMP
  inputBinding:
    prefix: --workingdir

outputs:
- id: out
  label: out
  type: File
  outputBinding:
    glob: generated.svs.vcf
    loadContents: false
- id: assembly
  label: assembly
  type: File
  outputBinding:
    glob: generated.assembled.bam
    loadContents: false
stdout: _stdout
stderr: _stderr

baseCommand: /opt/gridss/gridss.sh
arguments: []
id: gridss
