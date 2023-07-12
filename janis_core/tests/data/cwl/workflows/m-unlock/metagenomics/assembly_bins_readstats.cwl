#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

label: "Bin read mapping stats"

doc: |
    Table of general read mapping statistics of the bins and assembly
    
    ID
    Reads
    Assembly size
    Contigs
    n50
    Largest contig
    Mapped reads
    Bins
    Total bin size
    Binned
    Reads mapped to bins

requirements:
  - class: InlineJavascriptRequirement

hints:
  SoftwareRequirement:
    packages:
      python3:
        version: ["3.10.6"]
        specs: ["https://anaconda.org/conda-forge/python"]
      pysam:
        version: ["0.20.0"]
        specs: ["https://anaconda.org/bioconda/pysam"]       
  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/scripts:1.0.3

baseCommand: ["python3", "/scripts/metagenomics/assembly_bins_readstats.py"]

inputs:
  identifier:
    type: string
    doc: Identifier for this dataset used in this workflow
    label: identifier used
    inputBinding:
      prefix: --identifier
  
  binContigs:
    type: File
    doc: File containing bins names and their respective (assembly) contigs. Format contig<tab>bin_name
    label: binContigs file
    inputBinding:
      prefix: --binContigs
  bam_file:
    type: File
    doc: BAM file with reads mapped to the assembly
    label: BAM file
    inputBinding:
      prefix: --bam
  assembly:
    type: File
    doc: Assembly in fasta format
    label: Assembly
    inputBinding:
      prefix: --assembly

stdout: $(inputs.identifier)_binReadStats.tsv

outputs:
  binReadStats:
    type: File
    outputBinding:
      glob: $(inputs.identifier)_binReadStats.tsv

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
s:dateCreated: "2022-12-00"
s:license: https://spdx.org/licenses/Apache-2.0 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"

$namespaces:
  s: https://schema.org/
