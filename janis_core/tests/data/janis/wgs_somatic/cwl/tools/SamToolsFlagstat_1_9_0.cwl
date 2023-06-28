#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0
label: 'SamTools: Flagstat'
doc: |-
  Does a full pass through the input file to calculate and print statistics to stdout.

  Provides counts for each of 13 categories based primarily on bit flags in the FLAG field. Each category in the output is broken down into QC pass and QC fail. In the default output format, these are presented as "#PASS + #FAIL" followed by a description of the category.

  The first row of output gives the total number of reads that are QC pass and fail (according to flag bit 0x200). For example:

  122 + 28 in total (QC-passed reads + QC-failed reads)

  Which would indicate that there are a total of 150 reads in the input file, 122 of which are marked as QC pass and 28 of which are marked as "not passing quality controls"

  Following this, additional categories are given for reads which are:

  secondary     0x100 bit set

  supplementary     0x800 bit set

  duplicates     0x400 bit set

  mapped     0x4 bit not set

  paired in sequencing     0x1 bit set

  read1     both 0x1 and 0x40 bits set

  read2     both 0x1 and 0x80 bits set

  properly paired     both 0x1 and 0x2 bits set and 0x4 bit not set

  with itself and mate mapped     0x1 bit set and neither 0x4 nor 0x8 bits set

  singletons     both 0x1 and 0x8 bits set and bit 0x4 not set

  And finally, two rows are given that additionally filter on the reference name (RNAME), mate reference name (MRNM), and mapping quality (MAPQ) fields:

  with mate mapped to a different chr     0x1 bit set and neither 0x4 nor 0x8 bits set and MRNM not equal to RNAME

  with mate mapped to a different chr (mapQ>=5)     0x1 bit set and neither 0x4 nor 0x8 bits set and MRNM not equal to RNAME and MAPQ >= 5)

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: quay.io/biocontainers/samtools:1.9--h8571acd_11

inputs:
- id: bam
  label: bam
  type: File
  inputBinding:
    position: 10
- id: threads
  label: threads
  doc: Number of BAM compression threads to use in addition to main thread [0].
  type:
  - int
  - 'null'
  inputBinding:
    prefix: -@
    position: 5

outputs:
- id: out
  label: out
  type: stdout
stdout: _stdout
stderr: _stderr

baseCommand:
- samtools
- flagstat
arguments: []
id: SamToolsFlagstat
