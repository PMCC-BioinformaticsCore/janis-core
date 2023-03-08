cwlVersion: v1.0
class: CommandLineTool
id: gridss
label: gridss

doc: |
  Genomic Rearrangement IDentification Software Suite

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - $(inputs.input)
      - $(inputs.reference_genome)
      - $(inputs.blacklist)

hints:
  DockerRequirement:
    dockerPull: 'quay.io/biocontainers/gridss:2.9.3--0'

baseCommand: [gridss]

arguments:
  - position: 2
    prefix: '--output'
    valueFrom: $(inputs.input.nameroot.split('.')[0]).gridss.raw.vcf.gz
  - position: 5
    prefix: '--jar'
    valueFrom: "/usr/local/share/gridss-2.9.3-0/gridss.jar"

inputs:
  input:
    type: File
    secondaryFiles:
      - .bai
    inputBinding:
      position: 7
  reference_genome:
    type: File
    inputBinding:
      prefix: '--reference'
      position: 1
  assemblyFilename:
    type:
      - string?
      - "null"
    default: "gridss.assembly.bam"
    inputBinding:
      prefix: '--assembly'
      position: 3
  threads:
    type: int?
    default: 8
    inputBinding:
      prefix: '--threads'
      position: 4
  blacklist:
    type: File?
    inputBinding:
      prefix: '--blacklist'
      position: 6

outputs:
  output:
    type: File
    outputBinding:
      glob: "*.gridss.raw.vcf.gz"
      loadContents: false
