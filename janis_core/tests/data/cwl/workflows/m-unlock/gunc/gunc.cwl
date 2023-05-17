#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

label: GUNC

doc: the Genome UNClutterer, detection of chimerism and contamination in prokaryotic genomes.

requirements:
  - class: InlineJavascriptRequirement

hints:
  SoftwareRequirement:
    packages:
      gunc :
        version: ["1.0.5"]
        specs: ["https://anaconda.org/bioconda/gunc"]
  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/gunc:1.0.5

baseCommand: [ "gunc" ]

inputs:
  identifier:
    type: string
    doc: Identifier for this dataset used in this workflow
    label: identifier used
  threads:
    type: int?
    label: Number of threads to use
    default: 1
    inputBinding:
      position: 0
      prefix: -thread
  bin_dir:
    type: File
    doc: Input assembly in fasta format
    label: Input assembly
    inputBinding:
      prefix: --input_dir
  database:
    type: File
    doc: Abundances file
    label: Abundances
    inputBinding:
      prefix: --db_file

arguments:
  - prefix: "--out_dir"
    valueFrom: $(inputs.identifier)_GUNC

outputs:
# TODO

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
s:dateCreated: "2022-09-00"
s:license: https://spdx.org/licenses/Apache-2.0
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"

$namespaces:
  s: https://schema.org/
