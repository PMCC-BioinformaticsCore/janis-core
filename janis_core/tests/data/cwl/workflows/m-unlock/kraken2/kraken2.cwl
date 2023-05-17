#!/usr/bin/env cwltool

cwlVersion: v1.2
class: CommandLineTool

baseCommand: [ kraken2 ]

label: "Kraken2"
doc: |
    Kraken2 metagenomics taxomic read classification.
    
    Updated databases available at: https://benlangmead.github.io/aws-indexes/k2 (e.g. PlusPF-8)
    Original db: https://ccb.jhu.edu/software/kraken2/index.shtml?t=downloads

requirements:
  - class: InlineJavascriptRequirement

hints:
  SoftwareRequirement:
    packages:
      kraken2:
        version: ["2.1.2"]
        specs: ["https://anaconda.org/bioconda/kraken2"]

  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/kraken2:2.1.2

arguments:
  - valueFrom: $(inputs.identifier)_$(inputs.database.path.split( '/' ).pop())_kraken2.txt
    prefix: --output
  - valueFrom: $(inputs.identifier)_$(inputs.database.path.split( '/' ).pop())_kraken2_report.txt
    prefix: --report
  - "--report-zero-counts"
  - "--use-names"

inputs:
  threads:
    type: int?
    default: 1
    inputBinding:
      prefix: --threads
  identifier:
    type: string
    doc: Identifier for this dataset used in this workflow
    label: identifier used
  database:
    type: Directory
    label: Database
    doc: Database location of kraken2
    inputBinding:
      prefix: --db

  forward_reads:
    type: File?
    label: Forward reads
    doc: Illumina forward read file
    inputBinding:
      position: 100
  reverse_reads:
    type: File?
    label: Reverse reads
    doc: Illumina reverse read file
    inputBinding:
      position: 101
  paired_end:
    type: boolean?
    label: Paired end
    doc: Data is paired end (separate files)
    inputBinding:
      position: 2
      prefix: "--paired"
    default: false

  nanopore_reads:
    type: File?
    label: Nanopore reads
    doc: Oxford Nanopore Technologies reads in FASTQ
    inputBinding:
      position: 102

  confidence:
    type: float?
    label: Confidence
    doc: Confidence score threshold (default 0.0) must be in [0, 1]
    inputBinding:
      position: 4
      prefix: --confidence
  
  gzip:
    type: boolean?
    doc: "input data is gzip compressed"
    inputBinding:
      position: 3
      prefix: '--gzip-compressed'
    default: false
  bzip2:
    type: boolean
    doc: "input data is gzip compressed"
    inputBinding:
      position: 3
      prefix: '--bzip2-compressed'
    default: false

outputs:
  standard_report: 
    type: File
    outputBinding:
      glob: $(inputs.identifier)_$(inputs.database.path.split( '/' ).pop())_kraken2.txt
  sample_report:
    type: File
    outputBinding:
      glob: $(inputs.identifier)_$(inputs.database.path.split( '/' ).pop())_kraken2_report.txt

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-5516-8391
    s:email: mailto:german.royvalgarcia@wur.nl
    s:name: Germán Royval
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
s:dateCreated: "2021-11-25"
s:dateModified: "2021-11-04"
s:license: https://spdx.org/licenses/Apache-2.0 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"

$namespaces:
  s: https://schema.org/