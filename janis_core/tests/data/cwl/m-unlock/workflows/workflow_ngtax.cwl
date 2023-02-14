#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: Workflow
requirements:
   - class: StepInputExpressionRequirement
   - class: InlineJavascriptRequirement
   - class: MultipleInputFeatureRequirement

label: Quality assessment, amplicon classification
doc: | 
  Workflow for quality assessment of paired reads and classification using NGTax 2.0. 
  In addition files are exported to their respective subfolders for easier data management in a later stage.
  Steps:  
      - FastQC (read quality control)
      - NGTax 2.0
      - Export module

inputs:
  forward_reads:
    type: File
    doc: forward sequence file locally
    label: forward reads
  reverse_reads:
    type: File?
    doc: reverse sequence file locally
    label: reverse reads
  forward_primer: 
    type: string
    doc: Forward primer used
    label: Forward primer
  reverse_primer:
    type: string?
    doc: Reverse primer used
    label: Reverse primer
  reference_db:
    type: File?
    doc: Reference database used in FASTA format
    label: Reference database
  rev_read_len: 
    type: int?
    doc: Read length of the reverse read
    label: Reverse read length
  for_read_len: 
    type: int
    doc: Read length of the reverse read
    label: Reverse read length
  sample:
    type: string
    doc: Name of the sample being analysed
    label: Sample name
  fragment:
    type: string
    doc: Subfragment that is being analysed (e.g. V1-V3 or V5-region)
    label: Subfragment name
  primersRemoved:
    type: boolean?
    doc: Wether the primers are removed or not from the input files
    label: Primers are removed
  metadata:
    type: File?
    doc: UNLOCK assay metadata file
    label: Metadata file

  destination:
    type: string?
    label: Output Destination
    doc: Optional Output destination used for cwl-prov reporting.

steps:
############################
  fastqc:
    run: ../fastqc/fastqc.cwl
    in:
      fastq: 
        source: [forward_reads, reverse_reads]
        linkMerge: merge_flattened
        pickValue: all_non_null
    out: [html_files]
#############################
  reads_to_folder:
    run: ../expressions/files_to_folder.cwl
    in:
      files: 
        source: [forward_reads, reverse_reads]
        linkMerge: merge_flattened
        pickValue: all_non_null
      destination: 
        valueFrom: $("reads")
    out:
      [results]
############################
  ngtax:
    run: ../ngtax/ngtax.cwl
    in:
      forward_primer: forward_primer
      reverse_primer: reverse_primer
      reference_db: reference_db
      folder: reads_to_folder/results
      rev_read_len: rev_read_len
      for_read_len: for_read_len
      sample: sample
      fragment: fragment
      primersRemoved: primersRemoved
    out: [biom, turtle]
#############################
  ngtax_to_tsv-fasta:
    run: ../ngtax/ngtax_to_tsv-fasta.cwl
    in:
        input: ngtax/turtle
        identifier: sample
        fragment: fragment
        metadata: metadata
    out:
      [picrust_fasta, picrust_tsv, physeq_asv, physeq_seq, physeq_tax, physeq_met]
############################
  fastqc_files_to_folder:
    run: ../expressions/files_to_folder.cwl
    in:
      files: 
        source: [fastqc/html_files]
      destination: 
        valueFrom: $("1_QualityControl")
    out:
      [results]
############################
  ngtax_files_to_folder:
    run: ../expressions/files_to_folder.cwl
    in:
      files:
        source: [ngtax/biom, ngtax/turtle]
      destination: 
        valueFrom: $("2_Classification")
    out:
      [results]
############################
  phyloseq_files_to_folder:
    run: ../expressions/files_to_folder.cwl
    in:
      files: 
        source: [ngtax_to_tsv-fasta/physeq_asv, ngtax_to_tsv-fasta/physeq_seq, ngtax_to_tsv-fasta/physeq_tax, ngtax_to_tsv-fasta/physeq_met]
        linkMerge: merge_flattened
      destination:
        valueFrom: $("3_PHYLOSEQ")
    out:
      [results]
############################

outputs:
  turtle:
    type: File
    doc: Used for other workflows
    outputSource: ngtax/turtle
  files_to_folder_fastqc:
    type: Directory
    outputSource: fastqc_files_to_folder/results
  files_to_folder_ngtax:
    type: Directory
    outputSource: ngtax_files_to_folder/results
  files_to_folder_phyloseq:
    type: Directory
    outputSource: phyloseq_files_to_folder/results
    

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
s:dateCreated: "2021-01-01"
s:dateModified: "2023-01-18"
s:license: https://spdx.org/licenses/Apache-2.0 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"

$namespaces:
  s: https://schema.org/