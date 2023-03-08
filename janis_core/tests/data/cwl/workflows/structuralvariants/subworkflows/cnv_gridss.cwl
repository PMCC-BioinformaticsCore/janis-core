cwlVersion: v1.0
class: Workflow
label: cnv_gridss

doc: |
  CNV GRIDSS calling

requirements:
  ScatterFeatureRequirement: {}

inputs:
  bams: {type: 'File[]'}
  samples: {type: File}
  reference_genome: {type: File}
  blacklist: {type: 'File?'}
  threads: {type: 'int?'}
  min_len: {type: string}
  max_len: {type: string}
  min_q: {type: string}
  enable_gridss: {type: boolean}

outputs:
  output: {type: 'File?', outputSource: gridss_merge/output}

steps:
  gridss:
    run: ../tools/gridss.cwl
    in:
      input: bams
      reference_genome: reference_genome
      blacklist: blacklist
      threads: threads
    scatter: input
    out: [output]

  structural_variants:
    run: ../tools/structural_variants.cwl
    in:
      input:
        source: gridss/output
    scatter: input
    out: [output]

  gridss_filter:
    run: ../tools/gridss_filter.cwl
    in:
      input:
        source: structural_variants/output
      samples: samples
      min_len: min_len
      max_len: max_len
      min_q: min_q
    scatter: input
    out: [output]

  gridss_union:
    run: ../tools/bedops_union.cwl
    in:
      input:
        source: gridss_filter/output
    out: [output]

  gridss_merge:
    run: ../tools/collapse.cwl
    in:
      input:
        source: gridss_union/output
    out: [output]
