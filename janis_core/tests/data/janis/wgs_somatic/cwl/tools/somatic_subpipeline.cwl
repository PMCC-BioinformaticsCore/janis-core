#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.0

requirements:
- class: InlineJavascriptRequirement
- class: StepInputExpressionRequirement
- class: ScatterFeatureRequirement
- class: SubworkflowFeatureRequirement

inputs:
- id: reads
  type:
    type: array
    items:
      type: array
      items: File
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
- id: cutadapt_adapters
  type:
  - File
  - 'null'
- id: gatk_intervals
  type:
    type: array
    items: File
- id: snps_dbsnp
  type: File
  secondaryFiles:
  - .tbi
- id: snps_1000gp
  type: File
  secondaryFiles:
  - .tbi
- id: known_indels
  type: File
  secondaryFiles:
  - .tbi
- id: mills_indels
  type: File
  secondaryFiles:
  - .tbi
- id: align_and_sort_sortsam_tmpDir
  doc: Undocumented option
  type:
  - string
  - 'null'

outputs:
- id: out_bam
  type: File
  secondaryFiles:
  - .bai
  outputSource: merge_and_mark/out
- id: out_fastqc_reports
  type:
    type: array
    items:
      type: array
      items: File
  outputSource: fastqc/out
- id: out_performance_summary
  type: File
  outputSource: performance_summary/performanceSummaryOut

steps:
- id: fastqc
  label: FastQC
  in:
  - id: reads
    source: reads
  scatter:
  - reads
  run: fastqc_v0_11_8.cwl
  out:
  - id: out
  - id: datafile
- id: getfastqc_adapters
  label: Parse FastQC Adaptors
  in:
  - id: fastqc_datafiles
    source: fastqc/datafile
  - id: cutadapt_adaptors_lookup
    source: cutadapt_adapters
  scatter:
  - fastqc_datafiles
  run: ParseFastqcAdaptors_v0_1_0.cwl
  out:
  - id: adaptor_sequences
- id: align_and_sort
  label: Align and sort reads
  in:
  - id: sample_name
    source: sample_name
  - id: reference
    source: reference
  - id: fastq
    source: reads
  - id: cutadapt_adapter
    source: getfastqc_adapters/adaptor_sequences
  - id: cutadapt_removeMiddle3Adapter
    source: getfastqc_adapters/adaptor_sequences
  - id: sortsam_tmpDir
    source: align_and_sort_sortsam_tmpDir
  scatter:
  - fastq
  - cutadapt_adapter
  - cutadapt_removeMiddle3Adapter
  scatterMethod: dotproduct
  run: BwaAligner_1_0_0.cwl
  out:
  - id: out
- id: merge_and_mark
  label: Merge and Mark Duplicates
  in:
  - id: bams
    source: align_and_sort/out
  - id: sampleName
    source: sample_name
  run: mergeAndMarkBams_4_1_3.cwl
  out:
  - id: out
- id: calculate_performancesummary_genomefile
  label: Generate genome for BedtoolsCoverage
  in:
  - id: reference
    source: reference
  run: GenerateGenomeFileForBedtoolsCoverage_v0_1_0.cwl
  out:
  - id: out
- id: performance_summary
  label: Performance summary workflow (whole genome)
  in:
  - id: bam
    source: merge_and_mark/out
  - id: sample_name
    source: sample_name
  - id: genome_file
    source: calculate_performancesummary_genomefile/out
  run: PerformanceSummaryGenome_v0_1_0.cwl
  out:
  - id: performanceSummaryOut
id: somatic_subpipeline
