#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0
label: BGZip
doc: |-
  bgzip – Block compression/decompression utility

  Bgzip compresses files in a similar manner to, and compatible with, gzip(1). The file is compressed 
  into a series of small (less than 64K) 'BGZF' blocks. This allows indexes to be built against the 
  compressed file and used to retrieve portions of the data without having to decompress the entire file.

  If no files are specified on the command line, bgzip will compress (or decompress if the -d option is used) 
  standard input to standard output. If a file is specified, it will be compressed (or decompressed with -d). 
  If the -c option is used, the result will be written to standard output, otherwise when compressing bgzip 
  will write to a new file with a .gz suffix and remove the original. When decompressing the input file must 
  have a .gz suffix, which will be removed to make the output name. 
  Again after decompression completes the input file will be removed.

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: biodckrdev/htslib:1.2.1

inputs:
- id: file
  label: file
  doc: File to bgzip compress
  type: File
  inputBinding:
    position: 100
- id: outputFilename
  label: outputFilename
  type:
  - string
  - 'null'
  default: generated.gz
  inputBinding:
    position: 102
    valueFrom: $(inputs.file.basename).gz
    shellQuote: false
- id: offset
  label: offset
  doc: |-
    b: Decompress to standard output from virtual file position (0-based uncompressed offset). Implies -c and -d.
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --offset
- id: stdout
  label: stdout
  doc: 'c: Write to standard output, keep original files unchanged.'
  type: boolean
  default: true
  inputBinding:
    prefix: --stdout
- id: decompress
  label: decompress
  doc: 'd: Decompress.'
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --decompress
- id: force
  label: force
  doc: 'f: Overwrite files without asking.'
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --force
- id: help
  label: help
  doc: 'h: Displays a help message.'
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --help
- id: index
  label: index
  doc: |-
    i: Create a BGZF index while compressing. Unless the -I option is used, this will have the name of the compressed file with .gzi appended to it.
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --index
- id: indexName
  label: indexName
  doc: '-I: Index file name.'
  type:
  - File
  - 'null'
  inputBinding:
    prefix: --index-name
- id: compress
  label: compress
  doc: |-
    l: Compression level to use when compressing. From 0 to 9, or -1 for the default level set by the compression library. [-1]
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --compress
- id: reindex
  label: reindex
  doc: 'r: Rebuild the index on an existing compressed file.'
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --reindex
- id: rebgzip
  label: rebgzip
  doc: |-
    g: Try to use an existing index to create a compressed file with matching block offsets. Note that this assumes that the same compression library and level are in use as when making the original file. Don't use it unless you know what you're doing.
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --rebgzip
- id: size
  label: size
  doc: 's: Decompress INT bytes (uncompressed size) to standard output. Implies -c.'
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --size
- id: threads
  label: threads
  doc: '@: Number of threads to use [1].'
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --threads

outputs:
- id: out
  label: out
  type: File
  outputBinding:
    glob: $(inputs.file.basename).gz
    loadContents: false
stdout: _stdout
stderr: _stderr

baseCommand: bgzip
arguments:
- position: 101
  valueFrom: '>'
  shellQuote: false
id: bgzip
