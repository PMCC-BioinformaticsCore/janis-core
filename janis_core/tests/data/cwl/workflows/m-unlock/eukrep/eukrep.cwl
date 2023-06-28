#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

label: EukRep

doc: EukRep, Classification of Eukaryotic and Prokaryotic sequences from metagenomic datasets

requirements:
  - class: InlineJavascriptRequirement
  # - class: InitialWorkDirRequirement
    # listing:
    # - entry: "$({class: 'Directory', listing: []})"
     # entryname: "EukRep"
     # writable: true
    # - entryname: script.sh
      # entry: |-
          #!/bin/bash
          # source /root/miniconda/bin/activate
          # conda activate /unlock/infrastructure/conda/eukrep/eukrep_v0.6.7
          # EukRep $@

hints:
  SoftwareRequirement:
    packages:
      diamond :
        version: ["0.6.7"]
        specs: ["https://anaconda.org/bioconda/eukrep"]
  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/eukrep:0.6.7

# baseCommand: ["bash", "script.sh"] # see requirements
baseCommand: [ "EukRep" ]

inputs:
  identifier:
    type: string
    doc: Identifier for this dataset used in this workflow
    label: identifier used

  assembly:
    type: File
    doc: Input assembly in fasta format
    label: Input assembly
    inputBinding:
      prefix: -i

  stringency:
    type: string?
    doc: |
      {strict,balanced,lenient} Default is balanced.
      How stringent the algorithm is in identifying eukaryotic scaffolds. Strict has a lower false positive rate and true positive rate; vice verso for leneient.
    label: Algorithm stringency
    inputBinding:
      prefix: -m

  min_contig_size:
    type: string?
    label: Minumum contig length
    doc: Minimum sequence length cutoff for sequences to be included in prediction. Default is 3kb

arguments:
  - prefix: "-o"
    valueFrom: $(inputs.identifier)_EukRep.fasta

outputs:
  euk_fasta_out:
    type: File
    outputBinding:
      glob: $(inputs.identifier)_EukRep.fasta

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
s:dateCreated: "2022-06-00"
s:license: https://spdx.org/licenses/Apache-2.0
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"

$namespaces:
  s: https://schema.org/
