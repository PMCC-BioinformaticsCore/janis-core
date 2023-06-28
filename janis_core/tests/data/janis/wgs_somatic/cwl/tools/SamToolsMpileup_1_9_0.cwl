#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0
label: 'SamTools: Mpileup'
doc: |-
  Generate text pileup output for one or multiple BAM files. Each input file produces a separate group of pileup columns in the output.

  Samtools mpileup can still produce VCF and BCF output (with -g or -u), but this feature is deprecated and will be removed in a future release. Please use bcftools mpileup for this instead. (Documentation on the deprecated options has been removed from this manual page, but older versions are available online at <http://www.htslib.org/doc/>.)

  Note that there are two orthogonal ways to specify locations in the input file; via -r region and -l file. The former uses (and requires) an index to do random access while the latter streams through the file contents filtering out the specified regions, requiring no index. The two may be used in conjunction. For example a BED file containing locations of genes in chromosome 20 could be specified using -r 20 -l chr20.bed, meaning that the index is used to find chromosome 20 and then it is filtered for the regions listed in the bed file.

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: quay.io/biocontainers/samtools:1.9--h8571acd_11

inputs:
- id: illuminaEncoding
  label: illuminaEncoding
  doc: Assume the quality is in the Illumina 1.3+ encoding.
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --illumina1.3+
- id: countOrphans
  label: countOrphans
  doc: do not discard anomalous read pairs
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --count-orphans
- id: noBAQ
  label: noBAQ
  doc: disable BAQ (per-Base Alignment Quality)
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --no-BAQ
- id: adjustMQ
  label: adjustMQ
  doc: adjust mapping quality; recommended:50, disable:0 [0]
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --adjust-MQ
- id: maxDepth
  label: maxDepth
  doc: max per-file depth; avoids excessive memory usage [8000]
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --max-depth
- id: redoBAQ
  label: redoBAQ
  doc: recalculate BAQ on the fly, ignore existing BQs
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --redo-BAQ
- id: fastaRef
  label: fastaRef
  doc: ' skip unlisted positions (chr pos) or regions (BED)'
  type:
  - File
  - 'null'
  inputBinding:
    prefix: --fasta-ref
- id: excludeRG
  label: excludeRG
  doc: exclude read groups listed in FILE
  type:
  - File
  - 'null'
  inputBinding:
    prefix: --exclude-RG
- id: positions
  label: positions
  doc: skip unlisted positions (chr pos) or regions (BED)
  type:
  - File
  - 'null'
  inputBinding:
    prefix: --positions
- id: minBQ
  label: minBQ
  doc: Minimum base quality for a base to be considered [13]
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --min-BQ
- id: minMQ
  label: minMQ
  doc: skip alignments with mapQ smaller than INT [0]
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --min-MQ
- id: region
  label: region
  doc: region in which pileup is generated
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --region
- id: ignoreRG
  label: ignoreRG
  doc: ignore RG tags (one BAM = one sample)
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --ignore-RG
- id: inclFlags
  label: inclFlags
  doc: 'required flags: skip reads with mask bits unset []'
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --incl-flags
- id: exclFlags
  label: exclFlags
  doc: 'filter flags: skip reads with mask bits set [UNMAP,SECONDARY,QCFAIL,DUP]'
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --excl-flags
- id: ignoreOverlaps
  label: ignoreOverlaps
  doc: disable read-pair overlap detection
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --ignore-overlaps
- id: outputBP
  label: outputBP
  doc: output base positions on reads
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --output-BP
- id: outputMQ
  label: outputMQ
  doc: output mapping quality
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --output-MQ
- id: outputQNAME
  label: outputQNAME
  doc: output read names
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --output-QNAME
- id: allPositions
  label: allPositions
  doc: output all positions (including zero depth)
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: -a
- id: absolutelyAllPositions
  label: absolutelyAllPositions
  doc: output absolutely all positions, including unused ref. sequences
  type:
  - boolean
  - 'null'
- id: reference
  label: reference
  doc: Reference sequence FASTA FILE [null]
  type:
  - File
  - 'null'
  inputBinding:
    prefix: --reference
- id: bam
  label: bam
  type: File
  secondaryFiles:
  - .bai
  inputBinding:
    position: 10

outputs:
- id: out
  label: out
  type: stdout
stdout: _stdout
stderr: _stderr

baseCommand:
- samtools
- mpileup
arguments: []
id: SamToolsMpileup
