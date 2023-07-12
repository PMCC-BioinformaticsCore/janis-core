#!/usr/bin/env cwltool

cwlVersion: v1.2
class: CommandLineTool

hints:
  SoftwareRequirement:
    packages:
      rminionqc :
        version: ["1.4.2"]
        specs: ["https://anaconda.org/bioconda/r-minionqc"]

  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/r-minionqc:1.4.2


baseCommand: [ Rscript, /unlock/infrastructure/binaries/minionqc/MinIONQC.R ]

label: "MinIONQC for Quality Check of nanopore reads"
doc: |
    MinIONQC for Quality Check of nanopore reads

    Direct download: wget https://raw.githubusercontent.com/roblanf/minion_qc/master/MinIONQC.R -O MinIONQC.R

requirements:
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
    - entry: "$({class: 'Directory', listing: []})"
      entryname: "MinION_output"
      writable: true

arguments: 
  - valueFrom: "MinION_output"
    prefix: --outputdirectory

inputs:
  threads:
    type: int?
    default: 1
    inputBinding:
      prefix: --processors
  seq_summary:
    type: File
    doc: sequencing_summary.txt output from Guppy basecaller
    inputBinding:
      prefix: --input

outputs:
  qc_files: 
    type: File[]
    outputBinding:
      glob: MinION_output/*/*

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
s:license: https://spdx.org/licenses/Apache-2.0 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"

$namespaces:
  s: https://schema.org/
