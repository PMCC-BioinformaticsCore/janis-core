#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0
label: 'SamTools: View'
doc: |-
  Ensure SAMTOOLS.SORT is inheriting from parent metadata
          
  ---------------------------------------------------------------------------------------------------
      
  With no options or regions specified, prints all alignments in the specified input alignment file 
  (in SAM, BAM, or CRAM format) to standard output in SAM format (with no header).

  You may specify one or more space-separated region specifications after the input filename to 
  restrict output to only those alignments which overlap the specified region(s). 
  Use of region specifications requires a coordinate-sorted and indexed input file (in BAM or CRAM format).

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: quay.io/biocontainers/samtools:1.9--h8571acd_11

inputs:
- id: cramOutput
  label: cramOutput
  doc: Output in the CRAM format (requires -T).
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: -C
    position: 5
- id: compressedBam
  label: compressedBam
  doc: Enable fast BAM compression (implies -b).
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: '-1'
    position: 5
- id: uncompressedBam
  label: uncompressedBam
  doc: |-
    Output uncompressed BAM. This option saves time spent on compression/decompression and is thus preferred when the output is piped to another samtools command.
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: -u
    position: 5
- id: onlyOutputHeader
  label: onlyOutputHeader
  doc: Output the header only.
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: -H
    position: 5
- id: countAlignments
  label: countAlignments
  doc: |-
    Instead of printing the alignments, only count them and print the total number. All filter options, such as -f, -F, and -q, are taken into account.
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: -c
    position: 5
- id: writeAlignments
  label: writeAlignments
  doc: |-
    Write alignments that are not selected by the various filter options to FILE. When this option is used, all alignments (or all alignments intersecting the regions specified) are written to either the output file or this file, but never both.
  type:
  - File
  - 'null'
  inputBinding:
    prefix: -U
    position: 5
- id: inputTSV
  label: inputTSV
  doc: |-
    A tab-delimited FILE. Each line must contain the reference name in the first column and the length of the reference in the second column, with one line for each distinct reference. Any additional fields beyond the second column are ignored. This file also defines the order of the reference sequences in sorting. If you run: `samtools faidx <ref.fa>', the resulting index file <ref.fa>.fai can be used as this FILE.
  type:
  - File
  - 'null'
  inputBinding:
    prefix: -t
    position: 5
- id: onlyOverlapping
  label: onlyOverlapping
  doc: Only output alignments overlapping the input BED FILE [null].
  type:
  - File
  - 'null'
  inputBinding:
    prefix: -L
    position: 5
- id: useMultiRegionIterator
  label: useMultiRegionIterator
  doc: |-
    Use the multi-region iterator on the union of the BED file and command-line region arguments. This avoids re-reading the same regions of files so can sometimes be much faster. Note this also removes duplicate sequences. Without this a sequence that overlaps multiple regions specified on the command line will be reported multiple times.
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: -M
    position: 5
- id: outputAlignmentsInReadGroup
  label: outputAlignmentsInReadGroup
  doc: |-
    Output alignments in read group STR [null]. Note that records with no RG tag will also be output when using this option. This behaviour may change in a future release.
  type:
  - string
  - 'null'
  inputBinding:
    prefix: -r
    position: 5
- id: outputAlignmentsInFileReadGroups
  label: outputAlignmentsInFileReadGroups
  doc: |-
    Output alignments in read groups listed in FILE [null]. Note that records with no RG tag will also be output when using this option. This behaviour may change in a future release.
  type:
  - File
  - 'null'
  inputBinding:
    prefix: -R
    position: 5
- id: mapqThreshold
  label: mapqThreshold
  doc: Skip alignments with MAPQ smaller than INT [0].
  type:
  - int
  - 'null'
  inputBinding:
    prefix: -q
    position: 5
- id: outputAlignmentsInLibrary
  label: outputAlignmentsInLibrary
  doc: Only output alignments in library STR [null].
  type:
  - string
  - 'null'
  inputBinding:
    prefix: -l
    position: 5
- id: outputAlignmentsMeetingCIGARThreshold
  label: outputAlignmentsMeetingCIGARThreshold
  doc: |-
    Only output alignments with number of CIGAR bases consuming query sequence ≥ INT [0]
  type:
  - int
  - 'null'
  inputBinding:
    prefix: -m
    position: 5
- id: outputAlignmentsWithBitsSet
  label: outputAlignmentsWithBitsSet
  doc: |-
    Only output alignments with all bits set in INT present in the FLAG field. INT can be specified in hex by beginning with `0x' (i.e. /^0x[0-9A-F]+/) or in octal by beginning with `0' (i.e. /^0[0-7]+/) [0].
  type:
  - string
  - 'null'
  inputBinding:
    prefix: -f
    position: 5
- id: doNotOutputAlignmentsWithBitsSet
  label: doNotOutputAlignmentsWithBitsSet
  doc: |-
    Do not output alignments with any bits set in INT present in the FLAG field. INT can be specified in hex by beginning with `0x' (i.e. /^0x[0-9A-F]+/) or in octal by beginning with `0' (i.e. /^0[0-7]+/) [0].
  type:
  - string
  - 'null'
  inputBinding:
    prefix: -F
    position: 5
- id: doNotOutputAlignmentsWithAllBitsSet
  label: doNotOutputAlignmentsWithAllBitsSet
  doc: |-
    Do not output alignments with all bits set in INT present in the FLAG field. This is the opposite of -f such that -f12 -G12 is the same as no filtering at all. INT can be specified in hex by beginning with `0x' (i.e. /^0x[0-9A-F]+/) or in octal by beginning with `0' (i.e. /^0[0-7]+/) [0].
  type:
  - string
  - 'null'
  inputBinding:
    prefix: -G
    position: 5
- id: readTagToExclude
  label: readTagToExclude
  doc: Read tag to exclude from output (repeatable) [null]
  type:
  - string
  - 'null'
  inputBinding:
    prefix: -x
    position: 5
- id: collapseBackwardCIGAR
  label: collapseBackwardCIGAR
  doc: Collapse the backward CIGAR operation.
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: -B
    position: 5
- id: subsamplingProportion
  label: subsamplingProportion
  doc: |-
    Output only a proportion of the input alignments. This subsampling acts in the same way on all of the alignment records in the same template or read pair, so it never keeps a read but not its mate. The integer and fractional parts of the -s INT.FRAC option are used separately: the part after the decimal point sets the fraction of templates/pairs to be kept, while the integer part is used as a seed that influences which subset of reads is kept.
  type:
  - float
  - 'null'
  inputBinding:
    prefix: -s
    position: 5
- id: threads
  label: threads
  doc: Number of BAM compression threads to use in addition to main thread [0].
  type:
  - int
  - 'null'
  inputBinding:
    prefix: -@
    position: 5
- id: sam
  label: sam
  type: File
  inputBinding:
    position: 10
- id: reference
  label: reference
  doc: |-
    A FASTA format reference FILE, optionally compressed by bgzip and ideally indexed by samtools faidx. If an index is not present, one will be generated for you.
  type:
  - File
  - 'null'
  secondaryFiles:
  - .fai
  - .amb
  - .ann
  - .bwt
  - .pac
  - .sa
  - ^.dict
  inputBinding:
    prefix: -T
    position: 6
- id: outputFilename
  label: outputFilename
  doc: Output to FILE [stdout].
  type:
  - string
  - 'null'
  default: generated.bam
  inputBinding:
    prefix: -o
    position: 5
- id: regions
  label: regions
  doc: |-
    Region specifications after the input filename to restrict output to only those alignments which overlap the specified region(s). Use of region specifications requires a coordinate-sorted and indexed input file (in BAM or CRAM format)
  type:
  - type: array
    items: string
  - 'null'
  inputBinding:
    position: 11

outputs:
- id: out
  label: out
  type: File
  outputBinding:
    glob: generated.bam
    loadContents: false
stdout: _stdout
stderr: _stderr

baseCommand:
- samtools
- view
arguments:
- position: 2
  valueFrom: -S
- position: 3
  valueFrom: -h
- position: 4
  valueFrom: -b
id: SamToolsView
