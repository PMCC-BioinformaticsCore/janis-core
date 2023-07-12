cwlVersion: v1.0
class: CommandLineTool
id: fastqc
label: fastqc

doc: |
  Run fastqc on raw reads in FASTQ format (single or paired end)

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - $(inputs.fastq1)
      - $(inputs.fastq2)

hints:
  DockerRequirement:
    dockerPull: 'quay.io/biocontainers/fastqc:0.11.8--1'

baseCommand: [fastqc]

arguments:
  - valueFrom: $(runtime.outdir)
    prefix: '-o'
  - valueFrom: '--noextract'

inputs:
  fastq1:
    type: 'File[]'
    inputBinding:
      position: 2
  fastq2:
    type: 'File[]'
    inputBinding:
      position: 3
  threads:
    type: int?
    default: 2
    inputBinding:
      position: 1
      prefix: '--threads'

outputs:
  fastqc_zip:
    type:
      type: array
      items: File
    outputBinding:
      glob: '*_fastqc.zip'
  fastqc_html:
    type:
      type: array
      items: File
    outputBinding:
      glob: '*_fastqc.html'
