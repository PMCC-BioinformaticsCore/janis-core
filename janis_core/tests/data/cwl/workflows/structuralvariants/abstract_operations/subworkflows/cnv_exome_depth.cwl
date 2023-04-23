cwlVersion: v1.2
class: Workflow
label: cnv_exomedepth

doc: |
  CNV ExomeDepth calling

requirements:
  ScatterFeatureRequirement: {}

inputs:
  bams: {type: 'File[]'}
  samples: {type: File}
  reference_genome: {type: File}
  min_len: {type: string}
  max_len: {type: string}
  min_bf: {type: string}
  enable_exomeDepth: {type: boolean}

outputs:
  output: {type: 'File?', outputSource: exomedepth_merge/output}

steps:
  batch_parser:
    run: ../batch_parser.cwl
    in:
      input: bams
      samples: samples
    out: [output]

  exome_depth:
    run: ../exome_depth.cwl
    in:
      input: bams
      mapping:
        source: batch_parser/output
      reference_genome: reference_genome
    scatter: mapping
    out: [output]

  exomedepth_filter:
    run: ../exome_depth_filter.cwl
    in:
      input:
        source: exome_depth/output
      samples: samples
      min_len: min_len
      max_len: max_len
      min_bf: min_bf
    scatter: input
    out: [output]

  exomedepth_union:
    run: ../bedops_union.cwl
    in:
      input:
        source: exomedepth_filter/output
    out: [output]

  exomedepth_merge:
    run: ../collapse.cwl
    in:
      input:
        source: exomedepth_union/output
    out: [output]