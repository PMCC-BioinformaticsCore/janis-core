cwlVersion: v1.2
class: Workflow
label: indexing_bed

inputs:
  input: {type: File}

outputs:
  output: {type: File, outputSource: tabix/output}

steps:
  bgzip:
    run: ../bgzip.cwl
    in:
      input: input
    out: [output]

  tabix:
    run: ../tabix.cwl
    in:
      input:
        source: bgzip/output
    out: [output]
