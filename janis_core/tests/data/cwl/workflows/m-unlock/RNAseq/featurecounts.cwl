#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool

label: "Bowtie2 alignment"

doc: |
    Align reads to indexed genome. Stripped simple version; only paired end reads and sam output.

requirements:
 - class: InlineJavascriptRequirement

hints:
  SoftwareRequirement:
    packages:
      subread:
        version: ["2.0.1"]
        specs: ["https://anaconda.org/bioconda/subread"]
  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/subread:2.0.1

inputs:
  prefix:
    type: string?
    default: gene_counts_ftcounts

  threads:
    type: int?
    default: 1
    inputBinding:
      prefix: "-T"

  gtf:
    type: File
    inputBinding:
      prefix: "-a"  
  bam:
    type: File
    inputBinding:
      position: -1

arguments:
  - prefix: "-o"
    valueFrom: $(inputs.prefix)_FeatureCounts.txt

baseCommand: [ featureCounts ]

outputs:
  readcounts:
    type: File
    outputBinding:
      glob: "*FeatureCounts.txt"
  
  summary:
    type: File
    outputBinding:
      glob: "*.summary"

  overview:
    type: stderr
stderr: $(inputs.prefix)_FeatureCounts_overview.txt

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