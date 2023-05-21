#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

label: SemiBin

doc: Metagenomic binning with semi-supervised deep learning using information from reference genomes.

requirements:
  InlineJavascriptRequirement: {}
  NetworkAccess:
    networkAccess: true

hints:
  SoftwareRequirement:
    packages:
      semibin:
        version: ["1.4.0"]
        specs: ["https://anaconda.org/bioconda/semibin"]
  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/semibin:1.4.0

baseCommand: [ "SemiBin single_easy_bin" ]

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
      prefix: --threads
  assembly:
    type: File
    doc: Input assembly in fasta format
    label: Input assembly
    inputBinding:
      prefix: --input-fasta
  bam_file:
    type: File?
    doc: Mapped reads to assembly in sorted BAM format
    label: BAM file
    inputBinding:
      prefix: --input-bam
  metabat2_depth_file:
    type: File?
    doc: Contig depth file from MetaBAT2
    label: MetaBAT2 depths
    inputBinding:
      prefix: --depth-metabat2
  environment:
    type: string?
    doc: Built-in models (human_gut/dog_gut/ocean/soil/cat_gut/human_oral/mouse_gut/pig_gut/built_environment/wastewater/chicken_caecum/global)
    label: Environment
    inputBinding:
      prefix: --environment
  reference_database:
    type: Directory?
    doc: Reference Database data directory (usually, MMseqs2 GTDB)
    label: Reference Database
    inputBinding:
      prefix: --reference-db
  sequencing_type_longread:
    type: boolean?
    doc: An alternative binning algorithm for assemblies from long-read datasets. 
    label: Long read assembly
    inputBinding:
      prefix: --sequencing-type=long_read

arguments:
  - prefix: "-o"
    valueFrom: $(inputs.identifier)_SemiBin

outputs:
  recluster_bins:
    type: Directory
    label: Reclustered Bins
    doc: Directory of all reconstructed bins after reclustering
    outputBinding:
      glob: $(inputs.identifier)_SemiBin/output_recluster_bins
  output_bins:
    type: Directory
    label: Bins
    doc: Directory of all reconstructed bins before reclustering.
    outputBinding:
      glob: $(inputs.identifier)_SemiBin/output_prerecluster_bins
  info:
    type: File?
    label: Bins info
    doc: Info on (reclustered) bins (contig,nbs,n50 etc..)
    outputBinding:
      glob: $(inputs.identifier)_SemiBin/recluster_bins_info.tsv  
  mmseqs_contig_annotation:
    type: Directory?
    label: MMseqs annotation
    doc: MMseqs contig annotation
    outputBinding:
      glob: mmseqs_contig_annotation
  sample0:
    type: Directory?
    label: Markers
    doc: Directory with HMM marker hits
    outputBinding:
      glob: $(inputs.identifier)_SemiBin/sample0
  data_split:
    type: File?
    label: Training data
    doc: Data used in the training of deep learning model
    outputBinding:
      glob: $(inputs.identifier)_SemiBin/data_split.csv
  data:
    type: File
    label: Training data
    doc: Data used in the training of deep learning model
    outputBinding:
      glob: $(inputs.identifier)_SemiBin/data.csv
  model:
    type: File?
    label: Deep learning model
    doc: Saved semi-supervised deep learning model.
    outputBinding:
      glob: $(inputs.identifier)_SemiBin/model.h5
  coverage:
    type: File?
    label: Coverage data
    doc: Coverage data generated from depth file.
    outputBinding:
      glob: $(inputs.identifier)_SemiBin/"*_cov.csv"

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
s:dateModified: "2022-12-00"
s:license: https://spdx.org/licenses/Apache-2.0
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"

$namespaces:
  s: https://schema.org/
