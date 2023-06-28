cwlVersion: v1.0
class: Workflow
label: bwa_index

doc: |
  Modified from https://github.com/kids-first/kf-somatic-workflow/blob/master/sub_workflows/prepare_reference.cwl

inputs:
  generate_bwa_indexes: {type: 'boolean?'}
  reference_fasta: {type: File}
  reference_fai: {type: 'File?'}
  reference_amb: {type: 'File?'}
  reference_ann: {type: 'File?'}
  reference_bwt: {type: 'File?'}
  reference_pac: {type: 'File?'}
  reference_sa: {type: 'File?'}
  algoType: {type: 'string?'}

outputs:
  output: {type: File, secondaryFiles: ['.amb', '.ann', '.bwt', '.fai', '.pac', '.sa'], outputSource: bundle_secondaries/output}

steps:
  gunzip:
    run: ../tools/gunzip.cwl
    in:
      input_fasta: reference_fasta
    out: [fa]

  samtools_faidx:
    run: ../tools/samtools_faidx.cwl
    in:
      input_fasta:
        source: gunzip/fa
      input_index: reference_fai
    out: [fai]

  bwa_index:
    run: ../tools/bwa_index.cwl
    in:
      generate_bwa_indexes: generate_bwa_indexes
      input_fasta:
        source: gunzip/fa
      input_amb: reference_amb
      input_ann: reference_ann
      input_bwt: reference_bwt
      input_pac: reference_pac
      input_sa: reference_sa
      algoType: algoType
    out: [amb, ann, bwt, pac, sa]

  bundle_secondaries:
    run: ../tools/bundle_secondaryfiles.cwl
    in:
      primary_file:
        source: gunzip/fa
      secondary_files:
        source: [samtools_faidx/fai, bwa_index/amb, bwa_index/ann, bwa_index/bwt, bwa_index/pac, bwa_index/sa]
        linkMerge: merge_flattened
    out: [output]
