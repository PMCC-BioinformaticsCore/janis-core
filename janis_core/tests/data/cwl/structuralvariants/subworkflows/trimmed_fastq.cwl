cwlVersion: v1.0
class: Workflow
label: trimmed_fastq

doc: |
  Quality Control (raw data), Raw Data trimming and Quality Control (pre-processed)

requirements:
  ScatterFeatureRequirement: {}

inputs:
  fastq1: {type: 'File[]'}
  fastq2: {type: 'File[]'}
  threads_fastqc: {type: 'int?'}
  threads_fastp: {type: 'int?'}

outputs:
  fastqc_raw_zip: {type: 'File[]', outputSource: fastqc_raw/fastqc_zip}
  fastqc_raw_html: {type: 'File[]', outputSource: fastqc_raw/fastqc_html}
  paired_R1: {type: 'File[]', outputSource: fastp/paired_fastq1}
  paired_R2: {type: 'File[]', outputSource: fastp/paired_fastq2}
  unpaired_R1: {type: 'File[]', outputSource: fastp/unpaired_fastq1}
  unpaired_R2: {type: 'File[]', outputSource: fastp/unpaired_fastq2}
  html_report: {type: 'File[]', outputSource: fastp/html_report}
  json_report: {type: 'File[]', outputSource: fastp/json_report}
  fastqc_paired_zip: {type: 'File[]', outputSource: fastqc_pre_processed/fastqc_zip}
  fastqc_paired_html: {type: 'File[]', outputSource: fastqc_pre_processed/fastqc_html}

steps:
  fastqc_raw:
    run: ../tools/fastqc.cwl
    in:
      fastq1: fastq1
      fastq2: fastq2
      threads: threads_fastqc
    out: [fastqc_zip, fastqc_html]

  fastp:
    run: ../tools/fastp.cwl
    in:
      fastq1: fastq1
      fastq2: fastq2
      threads: threads_fastp
    scatter:
      - fastq1
      - fastq2
    scatterMethod: dotproduct
    out: [paired_fastq1, paired_fastq2, unpaired_fastq1, unpaired_fastq2, html_report, json_report]

  fastqc_pre_processed:
    run: ../tools/fastqc.cwl
    in:
      fastq1:
        source: fastp/paired_fastq1
      fastq2:
        source: fastp/paired_fastq2
      threads: threads_fastqc
    out: [fastqc_zip, fastqc_html]
