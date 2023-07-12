cwlVersion: v1.0
class: Workflow
label: cnv_manta

doc: |
  CNV Manta calling

requirements:
  ScatterFeatureRequirement: {}

inputs:
  bams: {type: 'File[]'}
  samples: {type: File}
  reference_genome: {type: File, secondaryFiles: ['.amb', '.ann', '.bwt', '.fai', '.pac', '.sa']}
  bed: {type: File, secondaryFiles: ['.tbi']}
  exome: {type: boolean}
  min_len: {type: string}
  max_len: {type: string}
  min_q: {type: string}
  enable_manta: {type: boolean}

outputs:
  output: {type: 'File?', outputSource: manta_merge/output}

steps:
  manta:
    run: ../tools/manta.cwl
    in:
      input: bams
      reference_genome: reference_genome
      regions: bed
      exome: exome
    scatter: input
    out: [output]

  svtools:
    run: ../tools/svtools.cwl
    in:
      input:
        source: manta/output
    scatter: input
    out: [output]

  manta_filter:
    run: ../tools/manta_filter.cwl
    in:
      input:
        source: svtools/output
      samples: samples
      min_len: min_len
      max_len: max_len
      min_q: min_q
    scatter: input
    out: [output]

  manta_union:
    run: ../tools/bedops_union.cwl
    in:
      input:
        source: manta_filter/output
    out: [output]

  manta_merge:
    run: ../tools/collapse.cwl
    in:
      input:
        source: manta_union/output
    out: [output]
