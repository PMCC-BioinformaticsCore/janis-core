#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

label: "CheckM"

doc: |
    CheckM provides a set of tools for assessing the quality of genomes recovered from isolates, single cells, or metagenomes

hints:
  SoftwareRequirement:
    packages:
      checkm-genome:
        version: ["1.2.1"]
        specs: ["https://anaconda.org/bioconda/checkm-genome"]

  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/checkm-genome:1.2.2

requirements:
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/checkm-genome:1.2.2
  InitialWorkDirRequirement:
    listing:
      - entry: "$({class: 'Directory', listing: []})"
        entryname: "/checkm_data"
        writable: true
        # Should be /venv/lib/python3.10/site-packages/checkm/DATA_CONFIG
      - entry: "$({class: 'Directory', listing: []})"
        entryname: "checkm_output"
        writable: true
      - entryname: script.sh
        entry: |-
          # !/bin/bash
          export CHECKM_DATA_PATH=/venv/checkm_data
          checkm lineage_wf $@

arguments:
  - position: 51
    valueFrom: "--reduced_tree"
  - position: 52
    prefix: "-f"
    valueFrom: $(inputs.identifier)_CheckM_report.txt
  - position: 54
    valueFrom: "$(inputs.identifier)_checkm"

baseCommand: ["bash","-x", "script.sh"]

inputs:
  identifier:
    type: string
    doc: Identifier for this dataset used in this workflow
    label: identifier used

  bin_dir:
    type: Directory
    doc: folder containing bins in fasta format from metagenomic binning
    label: bins folder
    inputBinding:
      position: 53

  threads:
    type: int?
    label: Number of threads to use
    default: 8
    inputBinding:
      position: 4
      prefix: -t

  fasta_extension:
    type: string?
    label: fasta file extension
    inputBinding:
      position: 5
      prefix: -x
    default: "fa"

outputs:
  checkm_out_folder:
    type: Directory
    outputBinding:
      glob: $(inputs.identifier)_checkm
  checkm_out_table:
    type: File
    outputBinding:
      glob: $(inputs.identifier)_CheckM_report.txt