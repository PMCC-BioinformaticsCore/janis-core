#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.0
label: Align and sort reads
doc: |-
  Align sorted bam with this subworkflow consisting of BWA Mem + SamTools + Gatk4SortSam

requirements:
- class: InlineJavascriptRequirement
- class: StepInputExpressionRequirement

inputs:
- id: sample_name
  type: string
- id: reference
  type: File
  secondaryFiles:
  - .fai
  - .amb
  - .ann
  - .bwt
  - .pac
  - .sa
  - ^.dict
- id: fastq
  type:
    type: array
    items: File
- id: cutadapt_adapter
  type:
  - type: array
    items: string
  - 'null'
- id: cutadapt_removeMiddle3Adapter
  type:
  - type: array
    items: string
  - 'null'
- id: cutadapt_front
  doc: |-
    (-g)  Sequence of an adapter ligated to the 5' end (paired data: of the first read). The adapter and any preceding bases are trimmed. Partial matches at the 5' end are allowed. If a '^' character is prepended ('anchoring'), the adapter is only found if it is a prefix of the read.
  type:
  - string
  - 'null'
- id: cutadapt_removeMiddle5Adapter
  doc: 5' adapter to be removed from second read in a pair.
  type:
  - string
  - 'null'
- id: cutadapt_qualityCutoff
  doc: |-
    (]3'CUTOFF, ]3'CUTOFF, -q)  Trim low-quality bases from 5' and/or 3' ends of each read before adapter removal. Applied to both reads if data is paired. If one value is given, only the 3' end is trimmed. If two comma-separated cutoffs are given, the 5' end is trimmed with the first cutoff, the 3' end with the second.
  type: int
  default: 15
- id: cutadapt_minimumLength
  doc: '(-m)  Discard reads shorter than LEN. Default: 0'
  type: int
  default: 50
- id: bwamem_markShorterSplits
  doc: Mark shorter split hits as secondary (for Picard compatibility).
  type: boolean
  default: true
- id: sortsam_sortOrder
  doc: |-
    The --SORT_ORDER argument is an enumerated type (SortOrder), which can have one of the following values: [unsorted, queryname, coordinate, duplicate, unknown]
  type: string
  default: coordinate
- id: sortsam_createIndex
  doc: Whether to create a BAM index when writing a coordinate-sorted BAM file.
  type: boolean
  default: true
- id: sortsam_validationStringency
  doc: |-
    Validation stringency for all SAM files read by this program. Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.The --VALIDATION_STRINGENCY argument is an enumerated type (ValidationStringency), which can have one of the following values: [STRICT, LENIENT, SILENT]
  type: string
  default: SILENT
- id: sortsam_maxRecordsInRam
  doc: |-
    When writing SAM files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort a SAM file, and increases the amount of RAM needed.
  type: int
  default: 5000000
- id: sortsam_tmpDir
  doc: Undocumented option
  type: string
  default: .

outputs:
- id: out
  type: File
  secondaryFiles:
  - .bai
  outputSource: sortsam/out

steps:
- id: cutadapt
  label: Cutadapt
  in:
  - id: fastq
    source: fastq
  - id: adapter
    source: cutadapt_adapter
  - id: front
    source: cutadapt_front
  - id: qualityCutoff
    source: cutadapt_qualityCutoff
  - id: minimumLength
    source: cutadapt_minimumLength
  - id: removeMiddle3Adapter
    source: cutadapt_removeMiddle3Adapter
  - id: removeMiddle5Adapter
    source: cutadapt_removeMiddle5Adapter
  run: cutadapt_2_1.cwl
  out:
  - id: out
- id: bwamem
  label: Bwa mem + Samtools View
  in:
  - id: reference
    source: reference
  - id: reads
    source: cutadapt/out
  - id: sampleName
    source: sample_name
  - id: markShorterSplits
    source: bwamem_markShorterSplits
  run: BwaMemSamtoolsView_0_7_17_1_9.cwl
  out:
  - id: out
- id: sortsam
  label: 'GATK4: SortSAM'
  in:
  - id: bam
    source: bwamem/out
  - id: sortOrder
    source: sortsam_sortOrder
  - id: createIndex
    source: sortsam_createIndex
  - id: maxRecordsInRam
    source: sortsam_maxRecordsInRam
  - id: tmpDir
    source: sortsam_tmpDir
  - id: validationStringency
    source: sortsam_validationStringency
  run: Gatk4SortSam_4_1_2_0.cwl
  out:
  - id: out
id: BwaAligner
