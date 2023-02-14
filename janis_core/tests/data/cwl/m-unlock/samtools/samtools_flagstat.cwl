#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

label: "samtools flagstat"
doc: "samtools flagstat - reports general alignment summary statistics"

requirements:
 - class: InlineJavascriptRequirement 
hints:
  SoftwareRequirement:
    packages:
      samtools:
        version: ["1.15.1"]
        specs: ["https://anaconda.org/bioconda/samtools"]

  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/samtools:1.15.1

inputs:
  identifier:
    type: string
    doc: Identifier for this dataset used in this workflow
    label: identifier used

  bam_file:
    type: File
    label: Bam file
    doc: (sorted) Bam file
    inputBinding:
      position: 1

  threads:
    type: int?
    label: Number of threads to use
    default: 2
    inputBinding:
      position: 0
      prefix: --threads

# baseCommand: ["/unlock/infrastructure/binaries/samtools/samtools_v1.15/bin/samtools","flagstat"]
baseCommand: [samtools, flagstat]

stdout: $(inputs.identifier)_flagstat.tsv

outputs:
  flagstat:
    type: File
    outputBinding:
      glob: $(inputs.identifier)_flagstat.tsv

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
s:dateCreated: "2021-00-00"
s:dateModified: "2022-02-22"
s:license: https://spdx.org/licenses/Apache-2.0 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"

$namespaces:
  s: https://schema.org/
