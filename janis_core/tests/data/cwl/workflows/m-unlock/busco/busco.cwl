#!/usr/bin/env cwltool

cwlVersion: v1.2
class: CommandLineTool

label: "BUSCO"
doc: |
    Based on evolutionarily-informed expectations of gene content of near-universal single-copy orthologs, 
    BUSCO metric is complementary to technical metrics like N50.

requirements:
  InlineJavascriptRequirement: {}
  NetworkAccess: 
    networkAccess: $(inputs.busco_data !== undefined)

hints:
  SoftwareRequirement:
    packages:
      busco:
        version: ["5.4.4"]
        specs: ["https://anaconda.org/bioconda/busco"]
  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/busco:5.4.4

baseCommand: ["busco"]

inputs:
  threads:
    type: int
    label: Number of threads
    default: 1
    inputBinding:
      prefix: --cpu

  # input:
  #     type:
  #       - type: record
  #         name: sequence_file
  #         fields:
  #           sequence_file:
  #             type: File
  #             label: Input fasta file
  #             doc: Input sequence file in FASTA format. Can be an assembled genome or transcriptome (DNA), or protein sequences from an annotated gene set. Also possible to use a path to a directory containing multiple input files.
  #             inputBinding:
  #               prefix: --in
  #       - type: record
  #         name: sequence_files
  #         fields:
  #           sequence_files:
  #             type: File
  #             label: Input fasta files
  #             doc: Multiple input sequence files in FASTA format. Can be an assembled genome or transcriptome (DNA), or protein sequences from an annotated gene set. Also possible to use a path to a directory containing multiple input files.
  #             inputBinding:
  #               prefix: --in
  #       - type: record
  #         name: sequence_folder
  #         fields:
  #           sequence_folder:
  #             type: Directory
  #             label: Input folder
  #             doc: Input folder with sequence files in FASTA format. Can be an assembled genome or transcriptome (DNA), or protein sequences from an annotated gene set. Also possible to use a path to a directory containing multiple input files.
  #             inputBinding:
  #               prefix: --in

  sequence_file:
    type: File?
    label: Input fasta file
    doc: Input sequence file in FASTA format. Can be an assembled genome or transcriptome (DNA), or protein sequences from an annotated gene set. Also possible to use a path to a directory containing multiple input files.
    inputBinding:
      prefix: --in

  sequence_folder:
    type: Directory?
    label: Input folder
    doc: Input folder with sequence files in FASTA format. Can be an assembled genome or transcriptome (DNA), or protein sequences from an annotated gene set. Also possible to use a path to a directory containing multiple input files.
    inputBinding:
      prefix: --in


  busco_data:
    type: Directory?
    label: Dataset location
    doc: This assumes --offline mode. Specify local filepath for finding BUSCO dataset downloads
    inputBinding:
      prefix: --download_path

  mode:
    type: string
    label: Input molecule type
    doc: |
      Specify which BUSCO analysis mode to run.
      There are three valid modes:
      - geno or genome, for genome assemblies (DNA)
      - tran or transcriptome, for transcriptome assemblies (DNA)
      - prot or proteins, for annotated gene sets (protein)
    inputBinding:
        prefix: --mode

  auto-lineage:
    type: boolean?
    label: Auto-lineage detection
    doc: Run auto-lineage to find optimum lineage path
    inputBinding:
      prefix: --auto-lineage 

  lineage:
    type: string?
    label: Lineage
    doc: Specify the name of the BUSCO lineage to be used.
    inputBinding:
        prefix: --lineage_dataset
  
  auto-lineage-prok:
    type: boolean?
    label: Prokaryote auto-lineage detection
    doc: Run auto-lineage just on non-eukaryote trees to find optimum lineage path.
    inputBinding:
      prefix: --auto-lineage-prok

  auto-lineage-euk:
    type: boolean?
    label: Eukaryote auto-lineage detection
    doc: Run auto-placement just on eukaryote tree to find optimum lineage path.
    inputBinding:
      prefix: --auto-lineage-euk

  identifier:
    type: string
    label: Name of the output file
    doc: Give your analysis run a recognisable short name. Output folders and files will be labelled with this name.
    inputBinding:
        prefix: --out
  
  tar_output:
    type: boolean
    label: Compress output
    doc: Compress some subdirectories with many files to save space
    default: true
    inputBinding:
        prefix: --tar

arguments:
 - |
    ${
      if (inputs.busco_data){
        return '--offline';
      } else {
        return null;
      }
    }

outputs:
  logs:
    label: BUSCO logs folder
    type: Directory?
    outputBinding:
      glob: $(inputs.identifier)/logs
  # run_folders:
  #   label: run_*_odb10 folders
  #   type: Directory[]?
  #   outputBinding:
  #     glob: $(inputs.identifier)/run_*
  # run_folders_batch:
  #   label: Genome run folders
  #   doc: Genome run folders when in batch mode
  #   type: Directory[]?
  #   outputBinding:
  #     glob: $(inputs.identifier)/*.*fa*
  short_summaries:
    label: BUSCO short summary files
    type: File[]?
    outputBinding:
      glob: $(inputs.identifier)/short_summary.*
  batch_summary:
    label: Batch summary
    type: File?
    doc: Summary file when input is multiple files
    outputBinding:
      glob: $(inputs.identifier)/batch_summary.txt

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
s:dateCreated: "2022-01-01"
s:dateModified: "2022-02-28"
s:license: https://spdx.org/licenses/Apache-2.0 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"

$namespaces:
  s: https://schema.org/
