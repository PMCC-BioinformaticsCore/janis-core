#!/usr/bin/env cwltool

cwlVersion: v1.2
class: CommandLineTool

hints:
  SoftwareRequirement:
    packages:
      flye:
        version: ["2.9.1"]
        specs: ["https://anaconda.org/bioconda/flye"]

  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/flye:2.9.1

baseCommand: [flye]

label: Flye
doc: Flye De novo assembler for single molecule sequencing reads, with a focus in Oxford Nanopore Technologies reads

inputs:
  output_folder_name:
    type: string?
    label: Output folder name
    doc: Output folder name
    inputBinding:
      prefix: --out-dir
    default: "flye_output"

  nano_raw: # FASTQ read, not FAST5
    type: File?
    label: ONT reads raw
    doc: ONT regular reads in FASTQ format, pre-Guppy5 (<20% error)
    inputBinding:
      prefix: --nano-raw
  nano_corrected:
    type: File?
    label: ONT corrected
    doc: ONT reads in FASTQ format that were corrected with other methods (<3% error)
    inputBinding:
      prefix: --nano-corr
  nano_high_quality:
    type: File?
    label: ONT high quality
    doc: ONT high-quality reads in FASTQ format, Guppy5 SUP or Q20 (<5% error)
    inputBinding:
      prefix: --nano-hq
  pacbio_raw:
    type: File?
    label: PacBio reads raw
    doc: PacBio regular CLR  reads in FASTQ format, (<20% error)
    inputBinding:
      prefix: --pacbio-raw
  pacbio_corrected:
    type: File?
    label: PacBio reads corrected
    doc: PacBio  reads in FASTQ format, that were corrected with other methods (<3% error)
    inputBinding:
      prefix: --pacbio-corr
  pacbio_hifi:
    type: File?
    label: PacBio HiFi reads
    doc: PacBio HiFi  reads in FASTQ format, (<1% error)
    inputBinding:
      prefix: --pacbio-hifi

  threads:
    type: int?
    default: 1
    inputBinding:
      prefix: --threads
  polishing_iterations:
    label: Flye will carry out polishing multiple times as determined here
    type: int?
    default: 1
    inputBinding:
      prefix: --iterations
  metagenome:
    type: boolean?
    label: Metagenome
    doc: Set to true if assembling a metagenome
    default: false
    inputBinding:
      prefix: --meta
  debug_mode:
    type: boolean?
    label: Debug mode
    doc: Set to true to display debug output while running
    default: false
    inputBinding:
      prefix: --debug
  genome_size:
    type: string?
    label: Genome size
    doc: Estimated genome size (for example, 5m or 2.6g)
    inputBinding:
      prefix: --genome-size


outputs:
  # flye_outdir:
  #   label: Directory containing all output produced by flye
  #   type: Directory
  #   outputBinding:
  #     glob: $(inputs.output_folder_name)
  00_assembly:
    type: Directory
    outputBinding:
      glob: $(inputs.output_folder_name)/00-assembly
  10_consensus:
    type: Directory
    outputBinding:
      glob: $(inputs.output_folder_name)/10-consensus
  20_repeat:
    type: Directory
    outputBinding:
      glob: $(inputs.output_folder_name)/20-repeat
  30_contigger:
    type: Directory
    outputBinding:
      glob: $(inputs.output_folder_name)/30-contigger
  40_polishing:
    type: Directory
    outputBinding:
      glob: $(inputs.output_folder_name)/40-polishing
  assembly:
    label: Polished assembly created by flye, main output for after polishing with next tool
    type: File
    outputBinding:
      glob: $(inputs.output_folder_name)/assembly.fasta
  assembly_info:
    type: File
    outputBinding:
      glob: $(inputs.output_folder_name)/assembly_info.txt
  flye_log:
    type: File
    outputBinding:
      glob: $(inputs.output_folder_name)/flye.log
  params:
    type: File
    outputBinding:
      glob: $(inputs.output_folder_name)/params.json

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
s:dateCreated: "2021-11-29"
s:dateModified: "2023-01-10"
s:license: https://spdx.org/licenses/Apache-2.0 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"

$namespaces:
  s: https://schema.org/
