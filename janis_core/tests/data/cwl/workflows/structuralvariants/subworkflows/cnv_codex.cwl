cwlVersion: v1.0
class: Workflow
label: cnv_codex

doc: |
  CNV CODEX calling

requirements:
  ScatterFeatureRequirement: {}

inputs:
  bams: {type: 'File[]'}
  samples: {type: File}
  bed: {type: File}
  chromosome: {type: string}
  min_len: {type: string}
  max_len: {type: string}
  min_lratio: {type: string}
  enable_codex: {type: boolean}

outputs:
  output: {type: 'File?', outputSource: codex_merge/output}

steps:
  batch_parser:
    run: ../tools/batch_parser.cwl
    in:
      input: bams
      samples: samples
    out: [output]

  codex:
    run: ../tools/codex.cwl
    in:
      input: bams
      mapping:
        source: batch_parser/output
      bed: bed
      chromosome: chromosome
    scatter: mapping
    out: [output]

  codex_filter:
    run: ../tools/codex_filter.cwl
    in:
      input:
        source: codex/output
      samples: samples
      min_len: min_len
      max_len: max_len
      min_lratio: min_lratio
    scatter: input
    out: [output]

  codex_union:
    run: ../tools/bedops_union.cwl
    in:
      input:
        source: codex_filter/output
    out: [output]

  codex_merge:
    run: ../tools/collapse.cwl
    in:
      input:
        source: codex_union/output
    out: [output]
