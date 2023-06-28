#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

label: "sam to PE fastq reads"
doc: |
  Extracts Paired End (PE) reads that are mapped OR unmapped the reference from a sam file.

requirements:
 - class: ShellCommandRequirement
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
  threads:
    type: int?
    default: 1
    inputBinding:
      position: 3
      prefix: --thread
  sam:
    type: File
    doc: SAM file with reads
    label: sam file
    inputBinding:
      position: 4
  unmapped:
    type: boolean?
    doc: Extract unmapped reads
    label: Unmapped Reads
    inputBinding:
      position: 1
      prefix: "-f 4"
  mapped:
    type: boolean?
    doc: Extract mapped reads
    label: Unmapped Reads
    inputBinding:
      position: 1
      prefix: "-F 4"
    default: true

outputs:
  out_forward_reads: 
    type: File
    outputBinding:
      glob: $(inputs.identifier)_filtered_1.fq.gz
  out_reverse_reads: 
    type: File
    outputBinding:
      glob: $(inputs.identifier)_filtered_2.fq.gz

# baseCommand: ["/unlock/infrastructure/binaries/samtools/samtools_v1.15/bin/samtools","fastq"]
baseCommand: ["samtools", "fastq"]

arguments:
  - prefix: "-o"
    position: 1
    valueFrom: $(inputs.identifier)_filtered.fq.gz

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
s:dateModified: "2022-04-00"
s:license: https://spdx.org/licenses/Apache-2.0 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"

$namespaces:
  s: https://schema.org/
