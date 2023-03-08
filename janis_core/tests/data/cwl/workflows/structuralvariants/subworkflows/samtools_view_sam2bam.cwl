cwlVersion: v1.0
class: Workflow
label: samtools_view_sam2bam

requirements:
  ScatterFeatureRequirement: {}

inputs:
  paired: {type: 'File[]'}
  unpaired_R1: {type: 'File[]'}
  unpaired_R2: {type: 'File[]'}
  threads: {type: 'int?'}

outputs:
  output_paired: {type: 'File[]', outputSource: samtools_view_sam2bam_paired/output}
  output_unpairedR1: {type: 'File[]', outputSource: samtools_view_sam2bam_unpaired_1/output}
  output_unpairedR2: {type: 'File[]', outputSource: samtools_view_sam2bam_unpaired_2/output}

steps:
  samtools_view_sam2bam_paired:
    run: ../tools/samtools_view_sam2bam.cwl
    in:
      input: paired
      threads: threads
    scatter: input
    out: [output]

  samtools_view_sam2bam_unpaired_1:
    run: ../tools/samtools_view_sam2bam.cwl
    in:
      input: unpaired_R1
      threads: threads
    scatter: input
    out: [output]

  samtools_view_sam2bam_unpaired_2:
    run: ../tools/samtools_view_sam2bam.cwl
    in:
      input: unpaired_R2
      threads: threads
    scatter: input
    out: [output]
