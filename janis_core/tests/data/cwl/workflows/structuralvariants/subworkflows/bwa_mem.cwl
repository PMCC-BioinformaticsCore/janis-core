cwlVersion: v1.0
class: Workflow
label: bwa_mem

requirements:
  MultipleInputFeatureRequirement: {}

inputs:
  paired_R1: {type: File}
  paired_R2: {type: File}
  unpaired_R1: {type: File}
  unpaired_R2: {type: File}
  reference_fasta: {type: File, secondaryFiles: ['.amb', '.ann', '.bwt', '.fai', '.pac', '.sa']}
  threads: {type: 'int?'}
  read_group: {type: 'string?'}

outputs:
  output_paired: {type: File, outputSource: bwa_mem_paired_reads/output}
  output_unpairedR1: {type: File, outputSource: bwa_mem_unpaired_reads_1/output}
  output_unpairedR2: {type: File, outputSource: bwa_mem_unpaired_reads_2/output}

steps:
  bwa_mem_paired_reads:
    run: ../tools/bwa_mem_paired.cwl
    in:
      reads:
        source: [paired_R1, paired_R2]
        linkMerge: merge_flattened
      reference_genome: reference_fasta
      threads: threads
      read_group: read_group
    out: [output]

  bwa_mem_unpaired_reads_1:
    run: ../tools/bwa_mem_unpaired.cwl
    in:
      reads: unpaired_R1
      reference_genome: reference_fasta
      threads: threads
      read_group: read_group
    out: [output]

  bwa_mem_unpaired_reads_2:
    run: ../tools/bwa_mem_unpaired.cwl
    in:
      reads: unpaired_R2
      reference_genome: reference_fasta
      threads: threads
      read_group: read_group
    out: [output]
