#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0
label: Tabix
doc: |-
  tabix – Generic indexer for TAB-delimited genome position files

  Tabix indexes a TAB-delimited genome position file in.tab.bgz and creates an index file (in.tab.bgz.tbi or 
  in.tab.bgz.csi) when region is absent from the command-line. The input data file must be position sorted 
  and compressed by bgzip which has a gzip(1) like interface.

  After indexing, tabix is able to quickly retrieve data lines overlapping regions specified in the format 
  "chr:beginPos-endPos". (Coordinates specified in this region format are 1-based and inclusive.)

  Fast data retrieval also works over network if URI is given as a file name and in this case the 
  index file will be downloaded if it is not present locally.

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
- class: InitialWorkDirRequirement
  listing:
  - entry: $(inputs.inp)
- class: DockerRequirement
  dockerPull: biodckrdev/htslib:1.2.1

inputs:
- id: inp
  label: inp
  doc: |-
    File from which to create the index. The input data file must be position sorted and compressed by bgzip which has a gzip(1) like interface.
  type: File
  inputBinding:
    position: 8
- id: preset
  label: preset
  doc: |-
    -p: Input format for indexing. Valid values are: gff, bed, sam, vcf. This option should not be applied together with any of -s, -b, -e, -c and -0; it is not used for data retrieval because this setting is stored in the index file. [gff]
  type: string
  default: vcf
  inputBinding:
    prefix: --preset
    position: 2
- id: zeroBased
  label: zeroBased
  doc: |-
    -0: Specify that the position in the data file is 0-based (e.g. UCSC files) rather than 1-based.
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --zero-based
    position: 1
- id: begin
  label: begin
  doc: '-b: Column of start chromosomal position. [4]'
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --begin
    position: 4
- id: comment
  label: comment
  doc: '-c: Skip lines started with character CHAR. [#]'
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --comment
    position: 7
- id: csi
  label: csi
  doc: '-C: Produce CSI format index instead of classical tabix or BAI style indices.'
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --csi
    position: 1
- id: end
  label: end
  doc: |-
    -e: Column of end chromosomal position. The end column can be the same as the start column. [5]
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --end
    position: 5
- id: force
  label: force
  doc: '-f: Force to overwrite the index file if it is present.'
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --force
    position: 1
- id: minShift
  label: minShift
  doc: '-m: set minimal interval size for CSI indices to 2^INT [14]'
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --min-shift
    position: 1
- id: sequence
  label: sequence
  doc: |-
    -s: Column of sequence name. Option -s, -b, -e, -S, -c and -0 are all stored in the index file and thus not used in data retrieval. [1]
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --sequence
    position: 3
- id: skipLines
  label: skipLines
  doc: '-S: Skip first INT lines in the data file. [0]'
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --skip-lines
    position: 6
- id: printHeader
  label: printHeader
  doc: '-h: Print also the header/meta lines.'
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --print-header
    position: 1
- id: onlyHeader
  label: onlyHeader
  doc: '-H: Print only the header/meta lines.'
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --only-header
    position: 1
- id: listChroms
  label: listChroms
  doc: '-l: List the sequence names stored in the index file.'
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --list-chroms
    position: 1
- id: reheader
  label: reheader
  doc: '-r: Replace the header with the content of FILE'
  type:
  - File
  - 'null'
  inputBinding:
    prefix: --reheader
    position: 1
- id: regions
  label: regions
  doc: |-
    -R: Restrict to regions listed in the FILE. The FILE can be BED file (requires .bed, .bed.gz, .bed.bgz file name extension) or a TAB-delimited file with CHROM, POS, and, optionally, POS_TO columns, where positions are 1-based and inclusive. When this option is in use, the input file may not be sorted.
  type:
  - File
  - 'null'
  inputBinding:
    prefix: --regions
    position: 11
- id: targets
  label: targets
  doc: |-
    -T: Similar to -R but the entire input will be read sequentially and regions not listed in FILE will be skipped
  type:
  - File
  - 'null'
  inputBinding:
    prefix: --targets
    position: 11

outputs:
- id: out
  label: out
  type: File
  secondaryFiles:
  - .tbi
  outputBinding:
    glob: $(inputs.inp.basename)
    loadContents: false
stdout: _stdout
stderr: _stderr

baseCommand: tabix
arguments: []
id: tabix
