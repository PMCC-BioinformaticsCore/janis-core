cwlVersion: v1.2
class: Operation
id: fastp
label: fastp

requirements:
  DockerRequirement:
    dockerPull: 'quay.io/biocontainers/fastp:0.20.0--hdbcaa40_0'

inputs:
  fastq1:
    type: File
  fastq2:
    type: File
  threads:
    type: 'int?'
  cut_right:
    type: 'boolean?'
  cut_right_window_size:
    type: 'int?'
  cut_right_mean_quality:
    type: 'int?'
  trim_tail1:
    type: 'int?'
  length_required:
    type: 'int?'

outputs:
  paired_fastq1:
    type: File
  paired_fastq2:
    type: File
  unpaired_fastq1:
    type: File
  unpaired_fastq2:
    type: File
  html_report:
    type: File
  json_report:
    type: File