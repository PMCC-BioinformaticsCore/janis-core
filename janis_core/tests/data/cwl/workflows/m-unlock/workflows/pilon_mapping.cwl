#!/usr/bin/env cwltool
cwlVersion: v1.2
class: Workflow
requirements:
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement: {}
  MultipleInputFeatureRequirement: {}
  SubworkflowFeatureRequirement: {}

label: Metagenomics workflow
doc: |
    Workflow pilon assembly polishing
    Steps:
      - BBmap (Read mapping to assembly)
      - Pilon

outputs:
  pilon_polished_assembly:
    label: Polished genome
    type: File
    outputSource: pilon/pilon_polished_assembly
  vcf:
    label: VCF file
    doc: Compressed VCF file containing the changes.
    type: File
    outputSource: vcf_compress/outfile
  log:
    label: Pilon log
    doc: Pilon log
    type: File
    outputSource: pilon/pilon_log

inputs:
  identifier:
    type: string
    doc: Identifier for this dataset used in this workflow
    label: identifier used
  illumina_forward_reads:
    type: File
    doc: forward sequence file locally
    label: forward reads
  illumina_reverse_reads:
    type: File
    doc: reverse sequence file locally
    label: reverse reads
  assembly:
    type: File
    doc: Assembly in fasta format
    label: Assembly
  fixlist:
    type: string?
    label: Pilon fix list
    doc: A comma-separated list of categories of issues to try to fix
    default: "all"
  threads:
    type: int?
    doc: number of threads to use for computational processes
    label: number of threads
    default: 2
  memory:
    type: int?
    doc: Maximum memory usage in megabytes
    label: Maximum memory in MB
    default: 40000

  destination:
    type: string?
    label: Output Destination
    doc: Optional Output destination used for cwl-prov reporting.

steps:
#############################################
#### BBmap read mapping (illumina reads) for polishing
  readmapping_pilon:
    label: BBMap read mapping
    doc: Illumina read mapping using BBmap on assembled contigs
    run: ../bbmap/bbmap.cwl
    in:
      threads: threads
      identifier: identifier
      reference: assembly
      forward_reads: illumina_forward_reads
      reverse_reads: illumina_reverse_reads
      memory: memory
    out: [sam, stats, covstats, log]
#############################################
#### Convert sam file to sorted bam
  sam_to_sorted_bam:
    label: sam conversion to sorted bam
    doc: Sam file conversion to a sorted bam file
    run: ../samtools/sam_to_sorted-bam.cwl
    in:
      identifier: identifier
      sam: readmapping_pilon/sam
      threads: threads
    out: [sortedbam]
#############################################
#### Index bam file
  bam_index:
    label: samtools index
    doc: Index file generation for sorted bam file
    run: ../samtools/samtools_index.cwl
    in:
      bam_file: sam_to_sorted_bam/sortedbam
      threads: threads
    out: [bam_index]
#############################################
#### Create hybrid bam bai file for Pilon cwl
  expressiontool_bam_index:
    label: CWL hybrid bam/bai file
    run: ../samtools/expression_bam_index.cwl
    in:
      bam_file:
        source: sam_to_sorted_bam/sortedbam
      bam_index: bam_index/bam_index
    out: [hybrid_bamindex]
#############################################
#### Pilon
  pilon:
    label: pilon
    doc: Pilon draft assembly polishing with the mapped reads
    run: ../pilon/pilon.cwl
    in:
      identifier: identifier
      threads: threads
      memory: memory
      assembly: assembly
      bam_file: expressiontool_bam_index/hybrid_bamindex
      fixlist: fixlist
    out: [pilon_polished_assembly, pilon_vcf, pilon_log]

#############################################
#### Gzip compression of files
  vcf_compress:
    run: ../bash/pigz.cwl
    in:
      inputfile: pilon/pilon_vcf
      threads: threads
    out: [outfile]

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
s:dateCreated: "2022-04-00"
s:license: https://spdx.org/licenses/Apache-2.0 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"

$namespaces:
  s: https://schema.org/