#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

label: "Bin summary"

doc: |
    Creates a summary table of the bins and their quality and taxonomy.

    Columns are:
    Bin
    Contigs
    Size
    Largest_contig
    N50
    GC
    avgDepth
    GTDB-Tk_taxonomy
    BUSCO_Taxonomy
    BUSCO_score
    CheckM_Completeness
    CheckM_Contamination
    CheckM_Strain-heterogeneity    

requirements:
  - class: InlineJavascriptRequirement

hints:
  SoftwareRequirement:
    packages:
      python3:
        version: ["3.10.6"]
        specs: ["https://anaconda.org/conda-forge/python"]
      pandas:
        version: ["1.5.0"]
        specs: ["https://anaconda.org/conda-forge/pandas"]
  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/scripts:1.0.3

baseCommand: ["python3", "/scripts/metagenomics/bins_summary.py"]

inputs:
  identifier:
    type: string
    doc: Identifier for this dataset used in this workflow
    label: identifier used  
  bin_dir:
    type: Directory
    doc: Directory containing bins in fasta format from metagenomic binning
    label: Bins directory
    inputBinding:
      prefix: "--bin_dir"
  bin_depths:
    type: File
    doc: MetaBAT2 aggregateDepths file
    label: bin depths
    inputBinding:
      prefix: "--bin_depths"
  busco_batch:
    type: File?
    doc: Directory containing BUSCO reports
    label: BUSCO folder
    inputBinding:
      prefix: "--busco_batch"
  checkm:
    type: File?
    doc: CheckM report file
    label: CheckM report
    inputBinding:
      prefix: "--checkm"
  gtdbtk:
    type: File?
    doc: CheckM report file
    label: CheckM report
    inputBinding:
      prefix: "--gtdbtk"

arguments:
  - prefix: "--output"
    valueFrom: $(inputs.identifier)

outputs:
  bins_summary_table:
    type: File
    outputBinding:
      glob: $(inputs.identifier)_binSummary.tsv
  bin_contigs:
    type: File
    outputBinding:
      glob: $(inputs.identifier)_binContigs.tsv

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
s:dateModified: "2022-12-00"
s:license: https://spdx.org/licenses/Apache-2.0 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"

$namespaces:
  s: https://schema.org/
