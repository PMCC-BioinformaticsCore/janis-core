#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

label: "NGTax demultixplexing"

doc: |
    Runs NGTAX amplicon demultiplexing procedure

requirements:
 - class: InlineJavascriptRequirement

hints:
  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/ngtax:2.2.6

inputs:
  forward_primer: # "[AG]GGATTAGATACCC"
    type: string
    doc: Forward primer used
    label: The forward primer used
    inputBinding:
      prefix: -for_p
  reverse_primer: # "CGAC[AG][AG]CCATGCA[ACGT]CACCT"
    type: string
    doc: Reverse primer used
    label: The reverse primer used
    inputBinding:
      prefix: -rev_p
  forward_reads:
    type: File
    doc: Forward library file
    label: The forward library used
  reverse_reads: 
    type: File
    doc: Reverse library file
    label: The reverse library used
  mapping_file:
    type: File
    doc: Mapping file containing barcode information
    label: The mapping file
    inputBinding:
      prefix: -mapFile

baseCommand: ["-demultiplex", "-output", "demultiplexed"]

arguments:
  - prefix: "-fastQ"
    valueFrom: $(inputs.forward_reads.path),$(inputs.reverse_reads.path)

outputs:
  output:
    type: File[]
    outputBinding:
      glob: "demultiplexed/*"

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0001-8172-8981
    s:email: mailto:jasper.koehorst@wur.nl
    s:name: Jasper Koehorst

s:citation: https://m-unlock.nl
s:codeRepository: https://gitlab.com/m-unlock/cwl
s:dateCreated: "2020-00-00"
s:license: https://spdx.org/licenses/Apache-2.0 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"

$namespaces:
  s: https://schema.org/
