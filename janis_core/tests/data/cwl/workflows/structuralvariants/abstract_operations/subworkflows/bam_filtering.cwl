cwlVersion: v1.2
class: Workflow
label: bam_filtering

doc: |
  BAM filtering

inputs:
  input: {type: File}
  min_mapping_quality: {type: int}
  bits_set: {type: int}
  threads: {type: 'int?'}

outputs:
  output: {type: File, outputSource: samtools_index/output}

steps:
  samtools_view:
    run: ../samtools_view.cwl
    in:
      input: input
      min_mapping_quality: min_mapping_quality
      bits_set: bits_set
      threads: threads
    out: [output]

  samtools_index:
    run: ../samtools_index.cwl
    in:
      input:
        source: samtools_view/output
    out: [output]
