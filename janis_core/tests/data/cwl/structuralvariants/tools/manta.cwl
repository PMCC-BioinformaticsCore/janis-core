cwlVersion: v1.0
class: CommandLineTool
id: manta
label: manta

doc: |
  Manta Structural Variant Caller

requirements:
  ShellCommandRequirement: {}
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - $(inputs.input)
      - $(inputs.reference_genome)
      - $(inputs.regions)

hints:
  DockerRequirement:
    dockerPull: 'quay.io/biocontainers/manta:1.6.0--py27_0'

arguments:
  - position: 0
    valueFrom: "/usr/local/share/manta-1.6.0-0/bin/configManta.py"
    shellQuote: false
  - position: 2
    valueFrom: $(";{runDir}/runWorkflow.py".replace(/\{runDir\}/g, inputs.runDir))
    shellQuote: false

inputs:
  input:
    type: File
    secondaryFiles:
      - .bai
    inputBinding:
      prefix: '--bam'
      position: 1
      shellQuote: false
  reference_genome:
    type: File
    secondaryFiles:
      - .fai
      - .amb
      - .ann
      - .bwt
      - .pac
      - .sa
    inputBinding:
      prefix: '--referenceFasta'
      position: 1
      shellQuote: false
  exome:
    type: boolean
    inputBinding:
      prefix: "--exome"
      position: 1
      shellQuote: false
  regions:
    type: File
    secondaryFiles:
      - .tbi
    inputBinding:
      prefix: '--callRegions'
      position: 1
      shellQuote: false
  runDir:
    type: string?
    default: "generated"
    inputBinding:
      prefix: '--runDir'
      position: 1
      shellQuote: false

outputs:
  python:
    type: File
    outputBinding:
      glob: $("{runDir}/runWorkflow.py".replace(/\{runDir\}/g, inputs.runDir))
  output:
    type: File
    outputBinding:
      glob: |-
        $("{runDir}/results/variants/diploidSV.vcf.gz".replace(/\{runDir\}/g, inputs.runDir))
      outputEval: |
        ${
          self[0].basename = inputs.input.nameroot.split('.', 2)[0] + ".manta.raw.vcf.gz";
          return self[0]
        }
