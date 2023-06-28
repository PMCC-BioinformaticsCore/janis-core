#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0
label: FastQC
doc: |-
  FastQC is a program designed to spot potential problems in high througput sequencing datasets. It runs a set of analyses on one or more raw sequence files in fastq or bam format and produces a report which summarises the results.
  FastQC will highlight any areas where this library looks unusual and where you should take a closer look. The program is not tied to any specific type of sequencing technique and can be used to look at libraries coming from a large number of different experiment types (Genomic Sequencing, ChIP-Seq, RNA-Seq, BS-Seq etc etc).

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: quay.io/biocontainers/fastqc:0.11.8--2

inputs:
- id: reads
  label: reads
  type:
    type: array
    items: File
  inputBinding:
    position: 5
- id: outdir
  label: outdir
  doc: |-
    (-o) Create all output files in the specified output directory. Please note that this directory must exist as the program will not create it.  If this option is not set then the output file for each sequence file is created in the same directory as the sequence file which was processed.
  type: string
  default: .
  inputBinding:
    prefix: --outdir
- id: casava
  label: casava
  doc: |-
    Files come from raw casava output. Files in the same sample group (differing only by the group number) will be analysed as a set rather than individually. Sequences with the filter flag set in the header will be excluded from the analysis. Files must have the same names given to them by casava (including being gzipped and ending with .gz) otherwise they won't be grouped together correctly.
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --casava
- id: nano
  label: nano
  doc: |-
    Files come from naopore sequences and are in fast5 format. In this mode you can pass in directories to process and the program will take in all fast5 files within those directories and produce a single output file from the sequences found in all files.
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --nano
- id: nofilter
  label: nofilter
  doc: |-
    If running with --casava then don't remove read flagged by casava as poor quality when performing the QC analysis.
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --nofilter
- id: extract
  label: extract
  doc: |-
    If set then the zipped output file will be uncompressed in the same directory after it has been created.  By default this option will be set if fastqc is run in non-interactive mode.
  type: boolean
  default: true
  inputBinding:
    prefix: --extract
- id: java
  label: java
  doc: |-
    (-j) Provides the full path to the java binary you want to use to launch fastqc. If not supplied then java is assumed to be in your path.
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --java
- id: noextract
  label: noextract
  doc: |-
    Do not uncompress the output file after creating it.  You should set this option if you donot wish to uncompress the output when running in non-interactive mode. 
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --noextract
- id: nogroup
  label: nogroup
  doc: |-
    Disable grouping of bases for reads >50bp. All reports will show data for every base in the read. WARNING: Using this option will cause fastqc to crash and burn if you use it on really long reads, and your plots may end up a ridiculous size. You have been warned! 
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --nogroup
- id: format
  label: format
  doc: |-
    (-f) Bypasses the normal sequence file format detection and forces the program to use the specified format.  Valid formats are bam,sam,bam_mapped,sam_mapped and fastq 
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --format
- id: threads
  label: threads
  doc: |-
    (-t) Specifies the number of files which can be processed simultaneously. Each thread will be allocated 250MB of memory so you shouldn't run more threads than your available memory will cope with, and not more than 6 threads on a 32 bit machine
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --threads
    valueFrom: |-
      $([inputs.runtime_cpu, 1, 1].filter(function (inner) { return inner != null })[0])
- id: contaminants
  label: contaminants
  doc: |-
    (-c) Specifies a non-default file which contains the list of contaminants to screen overrepresented sequences against. The file must contain sets of named contaminants in the form name[tab]sequence.  Lines prefixed with a hash will be ignored.
  type:
  - File
  - 'null'
  inputBinding:
    prefix: --contaminants
- id: adapters
  label: adapters
  doc: |-
    (-a) Specifies a non-default file which contains the list of adapter sequences which will be explicity searched against the library. The file must contain sets of named adapters in the form name[tab]sequence. Lines prefixed with a hash will be ignored.
  type:
  - File
  - 'null'
  inputBinding:
    prefix: --adapters
- id: limits
  label: limits
  doc: |-
    (-l) Specifies a non-default file which contains a set of criteria which will be used to determine the warn/error limits for the various modules.  This file can also be used to selectively  remove some modules from the output all together. The format needs to mirror the default limits.txt file found in the Configuration folder.
  type:
  - File
  - 'null'
  inputBinding:
    prefix: --limits
- id: kmers
  label: kmers
  doc: |-
    (-k) Specifies the length of Kmer to look for in the Kmer content module. Specified Kmer length must be between 2 and 10. Default length is 7 if not specified. 
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --kmers
- id: quiet
  label: quiet
  doc: (-q) Supress all progress messages on stdout and only report errors.
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --quiet
- id: dir
  label: dir
  doc: |-
    (-d) Selects a directory to be used for temporary files written when generating report images.Defaults to system temp directory if not specified.
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --dir

outputs:
- id: out
  label: out
  type:
    type: array
    items: File
  outputBinding:
    glob: '*.zip'
    loadContents: false
- id: datafile
  label: datafile
  type:
    type: array
    items: File
  outputBinding:
    glob: '*/fastqc_data.txt'
    loadContents: false
stdout: _stdout
stderr: _stderr

baseCommand: fastqc
arguments: []
id: fastqc
