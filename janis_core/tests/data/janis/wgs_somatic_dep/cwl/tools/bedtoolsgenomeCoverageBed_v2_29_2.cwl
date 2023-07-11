#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0
label: 'BEDTools: genomeCoverageBed'
doc: |-
  bedtools genomecov computes histograms (default), per-base reports (-d) and BEDGRAPH (-bg) summaries of feature coverage (e.g., aligned sequences) for a given genome. Note: 1. If using BED/GFF/VCF, the input (-i) file must be grouped by chromosome. A simple sort -k 1,1 in.bed > in.sorted.bed will suffice. Also, if using BED/GFF/VCF, one must provide a genome file via the -g argument. 2. If the input is in BAM (-ibam) format, the BAM file must be sorted by position. Using samtools sort aln.bam aln.sorted will suffice.

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0

inputs:
- id: depth
  label: depth
  doc: |-
    Report the depth at each genome position (with one-based coordinates). Default behavior is to report a histogram.
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: -d
- id: depthZero
  label: depthZero
  doc: |-
    Report the depth at each genome position (with zero-based coordinates). Reports only non-zero positions. Default behavior is to report a histogram.
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: -dz
- id: BedGraphFormat
  label: BedGraphFormat
  doc: |-
    Report depth in BedGraph format. For details, see: genome.ucsc.edu/goldenPath/help/bedgraph.html
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: -bg
- id: BedGraphFormata
  label: BedGraphFormata
  doc: |-
    Report depth in BedGraph format, as above (-bg). However with this option, regions with zero coverage are also reported. This allows one to quickly extract all regions of a genome with 0  coverage by applying: 'grep -w 0$' to the output.
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: -bga
- id: split
  label: split
  doc: |-
    Treat 'split' BAM or BED12 entries as distinct BED intervals when computing coverage. For BAM files, this uses the CIGAR 'N' and 'D' operations to infer the blocks for computing coverage. For BED12 files, this uses the BlockCount, BlockStarts, and BlockEnds fields (i.e., columns 10,11,12).
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: -split
- id: strand
  label: strand
  doc: |-
    (STRING): can be + or -. Calculate coverage of intervals from a specific strand. With BED files, requires at least 6 columns (strand is column 6).
  type:
  - string
  - 'null'
  inputBinding:
    prefix: -strand
- id: pairEnd
  label: pairEnd
  doc: Calculate coverage of pair-end fragments. Works for BAM files only
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: -pc
- id: fragmentSize
  label: fragmentSize
  doc: |-
    Force to use provided fragment size instead of read length. Works for BAM files only
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: -fs
- id: du
  label: du
  doc: |-
    Change strand af the mate read (so both reads from the same strand) useful for strand specific. Works for BAM files only
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: -du
- id: fivePos
  label: fivePos
  doc: Calculate coverage of 5' positions (instead of entire interval).
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: '-5'
- id: threePos
  label: threePos
  doc: Calculate coverage of 3' positions (instead of entire interval).
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: '-3'
- id: max
  label: max
  doc: |-
    Combine all positions with a depth >= max into a single bin in the histogram. Irrelevant for -d and -bedGraph
  type:
  - int
  - 'null'
  inputBinding:
    prefix: -max
- id: scale
  label: scale
  doc: |-
    Scale the coverage by a constant factor. Each coverage value is multiplied by this factor before being reported. Useful for normalizing coverage by, e.g., reads per million (RPM). Default is 1.0; i.e., unscaled.
  type:
  - float
  - 'null'
  inputBinding:
    prefix: -scale
- id: trackline
  label: trackline
  doc: |-
    Adds a UCSC/Genome-Browser track line definition in the first line of the output. - See here for more details about track line definition: http://genome.ucsc.edu/goldenPath/help/bedgraph.html - NOTE: When adding a trackline definition, the output BedGraph can be easily uploaded to the Genome Browser as a custom track, BUT CAN NOT be converted into a BigWig file (w/o removing the first line).
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: -trackline
- id: trackopts
  label: trackopts
  doc: |-
    Writes additional track line definition parameters in the first line. - Example: -trackopts 'name="My Track" visibility=2 color=255,30,30' Note the use of single-quotes if you have spaces in your parameters.
  type:
  - string
  - 'null'
  inputBinding:
    prefix: -trackopts
- id: inputBam
  label: inputBam
  doc: |-
    Input bam file. Note: BAM _must_ be sorted by position. A 'samtools sort <BAM>' should suffice.
  type:
  - File
  - 'null'
  inputBinding:
    prefix: -ibam
- id: inputBed
  label: inputBed
  doc: |-
    Input bed file. Must be grouped by chromosome. A simple 'sort -k 1,1 <BED> > <BED>.sorted' will suffice.
  type:
  - File
  - 'null'
  inputBinding:
    prefix: -iBed
- id: inputFile
  label: inputFile
  doc: Input file, can be gff/vcf.
  type:
  - File
  - 'null'
  inputBinding:
    prefix: -i
- id: genome
  label: genome
  doc: |-
    Genome file. The genome file should tab delimited and structured as follows: <chromName><TAB><chromSize>.
  type:
  - File
  - 'null'
  inputBinding:
    prefix: -g

outputs:
- id: out
  label: out
  type: stdout
stdout: _stdout
stderr: _stderr

baseCommand:
- genomeCoverageBed
arguments: []
id: bedtoolsgenomeCoverageBed
