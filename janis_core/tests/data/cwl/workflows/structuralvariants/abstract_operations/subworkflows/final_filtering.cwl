cwlVersion: v1.2
class: Workflow
label: final_filtering

doc: |
  Final filtering

inputs:
  manta_input: {type: 'File?'}
  gridss_input: {type: 'File?'}
  exomeDepth_input: {type: 'File?'}
  codex_input: {type: 'File?'}

outputs:
  output: {type: File, outputSource: merge/output}

steps:
  union:
    run: ../bedops_union.cwl
    in:
      input:
        source:
          - manta_input
          - gridss_input
          - exomeDepth_input
          - codex_input
        linkMerge: merge_flattened
        pickValue: all_non_null
    out: [output]

  merge:
    run: ../merge_all.cwl
    in:
      input:
        source: union/output
    out: [output]