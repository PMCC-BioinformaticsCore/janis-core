#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

label: "jgi_summarize_bam_contig_depths"

doc: |
    Summarize contig read depth from bam file for metabat2 binning.

requirements:
 - class: InlineJavascriptRequirement

hints:
  SoftwareRequirement:
    packages:
      metabat2:
        version: ["2.15"]
        specs: ["https://anaconda.org/bioconda/metabat2"]
  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/metabat2:2.15

inputs:
  identifier:
    type: string
    doc: Name of the output file
    label: output file name

  bamFile:
    type: File
    inputBinding:
      position: 2

# baseCommand: [/unlock/infrastructure/binaries/MetaBAT/metabat_v2.12.1/jgi_summarize_bam_contig_depths]
baseCommand: [jgi_summarize_bam_contig_depths]

arguments:
  - position: 1
    prefix: '--outputDepth'
    valueFrom:  $(inputs.identifier)_contigDepths.tsv


outputs:
  depths:
    type: File
    outputBinding:
      glob: $(inputs.identifier)_contigDepths.tsv