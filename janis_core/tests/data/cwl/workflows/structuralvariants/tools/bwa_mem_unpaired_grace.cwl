requirements:
  InitialWorkDirRequirement:
    listing:
      - $(inputs.reads)
      - $(inputs.reference_genome)
hints:
  DockerRequirement:
    dockerPull: 'quay.io/biocontainers/bwa:0.7.17--h84994c4_5'

baseCommand: [bwa, mem]

inputs:
  reads:
  reference_genome:
  threads:
  read_group:

outputs:
  stdout


