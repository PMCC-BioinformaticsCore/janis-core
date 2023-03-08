cwlVersion: v1.0
class: Workflow
label: samtools_sort

requirements:
  ScatterFeatureRequirement: {}

inputs:
  paired: {type: 'File[]'}
  unpaired_R1: {type: 'File[]'}
  unpaired_R2: {type: 'File[]'}
  threads: {type: 'int?'}

outputs:
  output_paired: {type: 'File[]', outputSource: samtools_sort_paired/output}
  output_unpairedR1: {type: 'File[]', outputSource: samtools_sort_unpaired_1/output}
  output_unpairedR2: {type: 'File[]', outputSource: samtools_sort_unpaired_2/output}

steps:
  samtools_sort_paired:
    run: ../tools/samtools_sort.cwl
    in:
      input: paired
      threads: threads
    scatter: input
    out: [output]

  samtools_sort_unpaired_1:
    run: ../tools/samtools_sort.cwl
    in:
      input: unpaired_R1
      threads: threads
    scatter: input
    out: [output]

  samtools_sort_unpaired_2:
    run: ../tools/samtools_sort.cwl
    in:
      input: unpaired_R2
      threads: threads
    scatter: input
    out: [output]
