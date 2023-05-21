#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.0
label: Merge and Mark Duplicates

requirements:
- class: InlineJavascriptRequirement
- class: StepInputExpressionRequirement
- class: MultipleInputFeatureRequirement

inputs:
- id: bams
  type:
    type: array
    items: File
  secondaryFiles:
  - .bai
- id: createIndex
  type: boolean
  default: true
- id: maxRecordsInRam
  type: int
  default: 5000000
- id: sampleName
  type:
  - string
  - 'null'
- id: mergeSamFiles_useThreading
  doc: |-
    Option to create a background thread to encode, compress and write to disk the output file. The threaded version uses about 20% more CPU and decreases runtime by ~20% when writing out a compressed BAM file.
  type: boolean
  default: true
- id: mergeSamFiles_validationStringency
  doc: |-
    Validation stringency for all SAM files read by this program. Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.The --VALIDATION_STRINGENCY argument is an enumerated type (ValidationStringency), which can have one of the following values: [STRICT, LENIENT, SILENT]
  type: string
  default: SILENT

outputs:
- id: out
  type: File
  secondaryFiles:
  - .bai
  outputSource: markDuplicates/out

steps:
- id: mergeSamFiles
  label: 'GATK4: Merge SAM Files'
  in:
  - id: bams
    source: bams
  - id: sampleName
    source: sampleName
  - id: useThreading
    source: mergeSamFiles_useThreading
  - id: createIndex
    source: createIndex
  - id: maxRecordsInRam
    source: maxRecordsInRam
  - id: validationStringency
    source: mergeSamFiles_validationStringency
  run: Gatk4MergeSamFiles_4_1_3_0.cwl
  out:
  - id: out
- id: markDuplicates
  label: 'GATK4: Mark Duplicates'
  in:
  - id: bam
    source:
    - mergeSamFiles/out
    linkMerge: merge_nested
  - id: createIndex
    source: createIndex
  - id: maxRecordsInRam
    source: maxRecordsInRam
  run: Gatk4MarkDuplicates_4_1_3_0.cwl
  out:
  - id: out
  - id: metrics
id: mergeAndMarkBams
