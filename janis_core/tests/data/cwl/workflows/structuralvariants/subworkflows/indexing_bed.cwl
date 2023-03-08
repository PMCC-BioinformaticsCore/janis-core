cwlVersion: v1.0
class: Workflow
label: indexing_bed

inputs:
  input: {type: File}

outputs:
  output: {type: File, outputSource: tabix/output}

steps:
  bgzip:
    run: ../tools/bgzip.cwl
    in:
      input: input
    out: [output]

  tabix:
    run: ../tools/tabix.cwl
    in:
      input:
        source: bgzip/output
    out: [output]
