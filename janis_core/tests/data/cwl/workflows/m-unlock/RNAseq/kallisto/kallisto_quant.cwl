#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool

label: "kallisto quantification"

doc: |-
    Pseudoalignment with the tool kallisto
    https://github.com/common-workflow-library/bio-cwl-tools/tree/release/Kallisto

requirements:
  InlineJavascriptRequirement: {}
  LoadListingRequirement:
    loadListing: shallow_listing

hints:
  SoftwareRequirement:
    packages:
      kallisto:
        version: ["0.48.0"]
        specs: ["https://anaconda.org/bioconda/kallisto"]
  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/kallisto:0.48.0

inputs:
  threads:
    type: int?
    default: 1
    inputBinding:
      position: 0
      prefix: "--threads="
      separate: false

  identifier:
    type: string
    doc: prefix of the filename outputs
    default: "sampleX"

  index:
    type: string
    doc: "kallisto index file"
    inputBinding:
      prefix: '--index='
      separate: false
      # # Return kallisto index .idx file in the given folder.
      # valueFrom: $(self.listing[0].path)

  forward_reads:
    type:
     - File
    inputBinding:
      position: 100

  reverse_reads:
    type:
     - "null"
     - File
    inputBinding:
      position: 101

  isSingle:
    type: boolean
    inputBinding:
      position: 2
      prefix: "--single"
    default: false

  #Optional Inputs

  isBias:
    type: boolean?
    inputBinding:
      prefix: "--bias"

  isFusion:
    type: boolean?
    inputBinding:
      prefix: "--fusion"

  isSingleOverhang:
    type: boolean?
    inputBinding:
      prefix: "--single-overhang"
  
  FragmentLength:
    type: double?
    inputBinding:
      separate: false
      prefix: "--fragment-length="
  
  StandardDeviation:
    type: double?
    inputBinding:
      prefix: "--sd"
  
  BootstrapSamples:
    type: int?
    inputBinding:
      separate: false
      prefix: "--bootstrap-samples="
    default: 100
  
  Seed:
    type: int?
    inputBinding:
      prefix: "--seed"

#Using record inputs to create mutually exclusive inputs
  Strand:
    type:
      - "null"
      - type: record
        name: forward
        fields:
          forward:
              type: boolean
              inputBinding:
                prefix: "--fr-stranded"

      - type: record
        name: reverse
        fields:
          reverse:
            type: boolean
            inputBinding:
              prefix: "--rf-stranded"

  PseudoBam:
    type: boolean?
    inputBinding:
      prefix: "--pseudobam"

# #Using record inputs to create dependent inputs

  GenomeBam:
    type:
      - "null"
      - type: record
        name: genome_bam
        fields:
          genomebam:
            type: boolean
            inputBinding:
              prefix: "--genomebam"

          gtf:
            type: File
            inputBinding:
              prefix: "--gtf"

          chromosomes:
            type: File
            inputBinding:
              prefix: "--chromosomes"

arguments:
  - "--plaintext"
  - prefix: "--output-dir="
    separate: false
    valueFrom: $(inputs.identifier)_kallisto

# baseCommand: [/unlock/infrastructure/binaries/kallisto/kallisto_v0.48.0/kallisto, quant]
baseCommand: [kallisto, quant]

stderr: $(inputs.identifier)_kallisto_summary.txt

outputs:
  abundance.h5:
    type: File?
    outputBinding:
      glob: $(inputs.identifier)_kallisto/abundance.h5
  abundance.tsv:
    type: File
    outputBinding:
      glob: $(inputs.identifier)_kallisto/abundance.tsv
  run_info:
    type: File
    outputBinding:
      glob: $(inputs.identifier)_kallisto/run_info.json

  summary:
    type: stderr

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0001-8172-8981
    s:email: mailto:jasper.koehorst@wur.nl
    s:name: Jasper Koehorst
  - class: s:Person
    s:identifier: https://orcid.org/0000-0001-9524-5964
    s:email: mailto:bart.nijsse@wur.nl
    s:name: Bart Nijsse

s:citation: https://m-unlock.nl
s:codeRepository: https://gitlab.com/m-unlock/cwl
s:dateCreated: "2020-00-00"
s:license: https://spdx.org/licenses/Apache-2.0 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"


$namespaces:
  s: https://schema.org/