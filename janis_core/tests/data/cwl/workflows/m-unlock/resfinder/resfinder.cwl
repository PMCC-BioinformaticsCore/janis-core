#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool

label: "ResFinder 4"

doc: |
    ResFinder identifies acquired antimicrobial resistance genes in total or partial sequenced isolates of bacteria.

# python3 /scratch2/bart/AR/resfinder/run_resfinder.py 
# --acquired 
# --inputfastq read_1.fastq.gz read_2.fastq.gz 
# --outputPath firmtest3 
# --species 'other' 
# --min_cov 0.6
# --threshold 0.9
# --db_path_res resfinder/db_resfinder 
# --kmaPath resfinder/cge/kma/kma
# --blastPath ...

requirements:
 - class: InlineJavascriptRequirement

hints:
  SoftwareRequirement:
    packages:
      minimap2:
        version: ["4.1.11"]
        specs: ["https://anaconda.org/bioconda/resfinder"]
  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/resfinder:4.1.11

inputs:
  identifier:
    type: string
    doc: Identifier used for the output folder and/or files
    label: output identifier

  forward_reads:
    doc: forward sequence fastq file locally
    type: File
  reverse_reads:
    doc: reverse sequence fastq file locally
    type: File

  species:
    type: string
    doc: Species in the sample
    label: species
    inputBinding:
      prefix: '--species'
  min_cov:
    type: string
    doc: Minimum (breadth-of) coverage of ResFinder
    label: min_cov
    inputBinding:
      prefix: '--min_cov'
  threshold:
    type: string
    doc: Threshold for identity of ResFinder
    label: threshold
    inputBinding:
      prefix: '--threshold'

  db_path_res:
    type: string
    doc: Path to the databases for ResFinder
    label: db_path_res
    inputBinding:
      prefix: '--db_path_res'
  kmaPath:
    type: string
    doc: Path to KMA executable
    label: KMA path
    inputBinding:
      prefix: '--kmaPath'
  blastPath:
    type: string
    doc: Path to blastn executable
    label: blastn path
    inputBinding:
      prefix: '--blastPath'


arguments:
  - prefix: "--inputfastq"
    valueFrom: "$(inputs.forward_reads) $(inputs.reverse_reads)"
  - prefix: "--outputPath"
    valueFrom: $(inputs.identifier)_ResFinder

baseCommand: cat

outputs:
  outputfolder: 
    type: Directory
    outputBinding:
      glob: $(inputs.identifier)_ResFinder

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
s:license: https://spdx.org/licenses/Apache-2.0 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"


$namespaces:
  s: https://schema.org/
