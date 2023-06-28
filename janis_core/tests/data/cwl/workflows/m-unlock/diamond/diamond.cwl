#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

doc: |
  Diamond workflow implementation

requirements:
 - class: InlineJavascriptRequirement

hints:
  SoftwareRequirement:
    packages:
      diamond :
        version: ["2.0.15"]
        specs: ["https://anaconda.org/bioconda/diamond"]
  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/diamond:2.0.15
    
inputs:
  database:
    type: string
    doc: Path to the DIAMOND database file. Since v2.0.8, a BLAST database can also be used here.
    label: Database file path
    inputBinding:
      prefix: --db
  forward_reads:
    type: File
    doc: forward sequence file locally
    label: forward reads
    # inputBinding:
      # prefix: --query
  reverse_reads:
    type: File
    doc: reverse sequence file locally
    label: reverse reads
  # query:
  #   type: File[]
  #   doc: Path to the query input file in FASTA or FASTQ format (may be gzip compressed). Two files that contain the same number of sequences may be supplied when running in blastx mode.
  #   label: Input file(s)
  #   inputBinding:
  #     prefix: --query
  align:
    type: string
    doc: Align amino acid (blastp) or DNA (blastx) sequences against a protein reference database
    label: Blast type
    default: blastx
  blocksize:
    type: int
    doc: Block size in billions of sequence letters to be processed at a time. This is the main parameter for controlling the program’s memory and disk space usage. Bigger numbers will increase the use of memory and temporary disk space, but also improve performance. The program can be expected to use roughly six times this number of memory (in GB).
    label: Block size
    inputBinding:
      prefix: --block-size
    default: 12
  maxtargetseq:
    type: int
    doc: maximum number of target sequences to report alignments for (default=25)
    label: Max target sequences
    inputBinding:
      prefix: --max-target-seqs
    default: 25
  indexchunks:
    type: int
    doc: The number of chunks for processing the seed index. This option can be additionally used to tune the performance. The default value is -c4, while setting this parameter to -c1 instead will improve the performance at the cost of increased memory use. Note that the very-sensitive and ultra-sensitive modes use -c1 by default.
    label: The number of chunks
    inputBinding:
      prefix: --index-chunks
    default: 1
  # taxonlist:
  #   type: File
  #   doc: Comma-separated list of NCBI taxonomic IDs to filter the database by.
  #   label: Taxon list file
  #   inputBinding:
  #     prefix: --taxonlist
  matrix:
    type: string?
    doc: Scoring matrix. The following matrices are supported, with the default being BLOSUM62.
    label: scoring matrix
    inputBinding:
      prefix: --matrix
    default: BLOSUM62
  output:
    type: string?
    doc: Path to the output file. If this parameter is omitted, the results will be written to the standard output and all other program output will be suppressed.
    label: output file
    inputBinding:
      prefix: --out
  outfmt:
    type: string
    doc: Format of the output file. See the diamond manual for accepted output formats
    label: Output format
    inputBinding:
      prefix: --outfmt
    default: "100" # Diamond format
  threads:
    type: int?
    default: 3
    inputBinding:
      prefix: --threads
  identifier:
    type: string
    doc: Identifier for this dataset used in this workflow
    label: identifier used

baseCommand: [/unlock/infrastructure/binaries/diamond/diamond_v2.0.8.146/diamond]

arguments:
  - "$(inputs.align)"
  - "--query"
  - "$(inputs.forward_reads.path)" 
  - "$(inputs.reverse_reads.path)"
  - "--salltitles"

outputs:
  output_diamond:
    type: File
    outputBinding:
      glob: $(inputs.output).daa

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
s:license: https://spdx.org/licenses/Apache-2.0 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"

$namespaces:
  s: https://schema.org/
