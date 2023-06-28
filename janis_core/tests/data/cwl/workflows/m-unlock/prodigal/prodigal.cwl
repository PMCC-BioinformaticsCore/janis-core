#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool

label: "Prodigal"
doc: "Prokaryotic gene prediction using Prodigal"

hints:
  SoftwareRequirement:
    packages:
      prodigal:
        version: ["2.6.3"]
        specs: ["https://anaconda.org/bioconda/prodigal"]

  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/prodigal:2.6.3

baseCommand: [ prodigal ]

arguments:
  # - valueFrom: "sco" # What is the sco format?
  #   prefix: "-f"
  - valueFrom: $(inputs.input_fasta.nameroot).prodigal
    prefix: "-o"
  - valueFrom: $(inputs.input_fasta.nameroot).prodigal.ffn
    prefix: "-d"
  - valueFrom: $(inputs.input_fasta.nameroot).prodigal.faa
    prefix: "-a"

outputs:
  predicted_proteins_out:
    type: File
    outputBinding:
      glob: $(inputs.input_fasta.nameroot).prodigal
  predicted_proteins_ffn:
    type: File
    outputBinding:
      glob: $(inputs.input_fasta.nameroot).prodigal.ffn
  predicted_proteins_faa:
    type: File
    outputBinding:
      glob: $(inputs.input_fasta.nameroot).prodigal.faa

inputs:
  input_fasta:
    type: File
    inputBinding:
      prefix: "-i"

  meta_mode:
    type: boolean?
    doc: Input is a meta-genome
    inputBinding:
      prefix: -p
      valueFrom: meta

  single_mode:
    type: boolean?
    doc: Input is an isolate genome
    inputBinding:
      prefix: -p
      valueFrom: single

  # mode:
  #   type:
  #     - type: record
  #       name: single
  #       fields:
  #         single:
  #           type: boolean?
  #           inputBinding:
  #             prefix: -p
  #             valueFrom: single
  #     - type: record
  #       name: meta
  #       fields:
  #         meta:
  #           type: boolean?
  #           inputBinding:
  #             prefix: -p
  #             valueFrom: meta

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0001-8172-8981
    s:email: mailto:jasper.koehorst@wur.nl
    s:name: Jasper Koehorst
  - class: s:Person
    s:identifier: https://orcid.org/0000-0001-9524-5964
    s:email: mailto:bart.nijsse@wur.nl
    s:name: Bart Nijsse
  - class: s:Person
    s:identifier: https://orcid.org/0000-0001-6867-2039
    s:name:  Ekaterina Sakharova

s:copyrightHolder': EMBL - European Bioinformatics Institute
s:license': "https://www.apache.org/licenses/LICENSE-2.0"

s:citation: https://m-unlock.nl
s:codeRepository: https://gitlab.com/m-unlock/cwl
s:dateCreated: "2022-06-00"
s:dateModified: "2022-08-00"
s:license: https://spdx.org/licenses/Apache-2.0 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"
s:copyrightNotice: " Copyright < 2022 EMBL - European Bioinformatics Institute
    This file has been modified by UNLOCK - Unlocking Microbial Potential
"

$namespaces:
  s: https://schema.org/