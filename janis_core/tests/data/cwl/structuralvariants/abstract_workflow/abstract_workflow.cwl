cwlVersion: v1.2
class: Workflow
label: CNV_pipeline

requirements:
  InlineJavascriptRequirement: {}
  ScatterFeatureRequirement: {}
  MultipleInputFeatureRequirement: {}
  SubworkflowFeatureRequirement: {}

$namespaces:
  s: https://schema.org/

inputs:
  fastq1: {type: 'File[]', doc: "FASTQ 1 list of files"}
  fastq2: {type: 'File[]', doc: "FASTQ 2 list of files"}
  generate_bwa_indexes: {type: 'boolean?', doc: "enable generation of reference genome indexes"}
  reference_fasta: {type: File, doc: "reference genome"}
  reference_fai: {type: 'File?', doc: "fai index for reference genome"}
  reference_amb: {type: 'File?', doc: "amb index for reference genome"}
  reference_ann: {type: 'File?', doc: "ann index for reference genome"}
  reference_bwt: {type: 'File?', doc: "bwt index for reference genome"}
  reference_pac: {type: 'File?', doc: "pac index reference genome"}
  reference_sa: {type: 'File?', doc: "sa index reference genome"}
  bed: {type: File, secondaryFiles: ['.tbi?'], doc: "capture kit"}
  samples: {type: File, doc: "TXT file mapping cases ID, samples ID and batch"}
  threads_fastqc: {type: 'int?', doc: "number of threads for fastqc"}
  threads_fastp: {type: 'int?', doc: "number of threads for fastp"}
  threads_bwa_mem: {type: 'int?', doc: "number of threads for bwa-mem"}
  threads_samtools: {type: 'int?', doc: "number of threads for samtools"}
  threads_gridss: {type: 'int?', doc: "number of threads for GRIDSS"}
  bwt_algorithm: {type: 'string?', doc: "BWT construction algorithm: bwtsw or is (Default: auto)"}
  read_group: {type: 'string?', doc: "FASTA read group (used to fix BAM files)"}
  min_mapping_quality: {type: int, doc: "skip alignments with MAPQ smaller than this minimum"}
  bits_set: {type: int, doc: "skip output alignments with this bits set"}
  manta_exome: {type: boolean, doc: "provide appropriate settings for WES"}
  manta_min_len: {type: string, doc: "minimum CNV lenght for MANTA"}
  manta_max_len: {type: string, doc: "maximum CNV lenght for MANTA"}
  manta_min_q: {type: string, doc: "minimum CNV quality for MANTA"}
  blacklist: {type: 'File?', doc: "BED file containing regions to ignore"}
  gridss_min_len: {type: string, doc: "minimum CNV lenght for GRIDSS"}
  gridss_max_len: {type: string, doc: "maximum CNV lenght for GRIDSS"}
  gridss_min_q: {type: string, doc: "minimum CNV quality for GRIDSS"}
  exomeDepth_min_len: {type: string, doc: "minimum CNV lenght for EXOME_DEPTH"}
  exomeDepth_max_len: {type: string, doc: "maximum CNV lenght for EXOME_DEPTH"}
  exomeDepth_min_bf: {type: string, doc: "minimum CNV Bayes factor for EXOME_DEPTH"}
  chromosome: {type: string, doc: "chromosome for targeted sequencing for CODEX"} # FIXME: change codex.cwl to run all the chromosomes
  codex_min_len: {type: string, doc: "minimum CNV lenght for CODEX"}
  codex_max_len: {type: string, doc: "maximum CNV lenght for CODEX"}
  codex_min_lratio: {type: string, doc: "minimum CNV lratio for CODEX"}
  enable_manta: {type: boolean, doc: "execute MANTA predictions"}
  enable_gridss: {type: boolean, doc: "execute GRIDSS predictions"}
  enable_exomeDepth: {type: boolean, doc: "execute EXOME_DEPTH predictions"}
  enable_codex: {type: boolean, doc: "execute CODEX predictions"}

outputs:
  fastqc_raw_zip: {type: 'File[]', outputSource: trimmed_fastq/fastqc_raw_zip}
  fastqc_raw_html: {type: 'File[]', outputSource: trimmed_fastq/fastqc_raw_html}
  html_report: {type: 'File[]', outputSource: trimmed_fastq/html_report}
  json_report: {type: 'File[]', outputSource: trimmed_fastq/json_report}
  fastqc_paired_zip: {type: 'File[]', outputSource: trimmed_fastq/fastqc_paired_zip}
  fastqc_paired_html: {type: 'File[]', outputSource: trimmed_fastq/fastqc_paired_html}
  output_bam_filtering: {type: 'File[]', outputSource: bam_filtering/output}
  output_manta: {type: 'File?', outputSource: cnv_manta/output}
  output_gridss: {type: 'File?', outputSource: cnv_gridss/output}
  output_exomedepth: {type: 'File?', outputSource: cnv_exomedepth/output}
  output_codex: {type: 'File?', outputSource: cnv_codex/output}
  output_all: {type: File, outputSource: final_filtering/output}

steps:
  bwa_index:
    run: ../abstract_operations/subworkflows/bwa_index.cwl
    in:
      generate_bwa_indexes: generate_bwa_indexes
      reference_fasta: reference_fasta
      reference_fai: reference_fai
      reference_amb: reference_amb
      reference_ann: reference_ann
      reference_bwt: reference_bwt
      reference_pac: reference_pac
      reference_sa: reference_sa
      algoType: bwt_algorithm
    out: [output]

  trimmed_fastq:
    run: ../abstract_operations/subworkflows/trimmed_fastq.cwl
    in:
      fastq1: fastq1
      fastq2: fastq2
      threads_fastqc: threads_fastqc
      threads_fastp: threads_fastp
    out: [fastqc_raw_zip, fastqc_raw_html, paired_R1, paired_R2, unpaired_R1, unpaired_R2, html_report, json_report, fastqc_paired_zip, fastqc_paired_html]

  bwa_mem:
    run: ../abstract_operations/subworkflows/bwa_mem.cwl
    in:
      paired_R1:
        source: trimmed_fastq/paired_R1
      paired_R2:
        source: trimmed_fastq/paired_R2
      unpaired_R1:
        source: trimmed_fastq/unpaired_R1
      unpaired_R2:
        source: trimmed_fastq/unpaired_R2
      reference_fasta:
        source: bwa_index/output
      threads: threads_bwa_mem
      read_group: read_group
    scatter:
      - paired_R1
      - paired_R2
      - unpaired_R1
      - unpaired_R2
    scatterMethod: dotproduct
    out: [output_paired, output_unpairedR1, output_unpairedR2]

  samtools_view_sam2bam:
    run: ../abstract_operations/subworkflows/samtools_view_sam2bam.cwl
    in:
      paired:
        source: bwa_mem/output_paired
      unpaired_R1:
        source: bwa_mem/output_unpairedR1
      unpaired_R2:
        source: bwa_mem/output_unpairedR2
      threads: threads_samtools
    out: [output_paired, output_unpairedR1, output_unpairedR2]

  samtools_sort:
    run: ../abstract_operations/subworkflows/samtools_sort.cwl
    in:
      paired:
        source: samtools_view_sam2bam/output_paired
      unpaired_R1:
        source: samtools_view_sam2bam/output_unpairedR1
      unpaired_R2:
        source: samtools_view_sam2bam/output_unpairedR2
      threads: threads_samtools
    out: [output_paired, output_unpairedR1, output_unpairedR2]

  samtools_merge:
    run: ../abstract_operations/samtools_merge.cwl
    in:
      paired:
        source: samtools_sort/output_paired
      unpaired_R1:
        source: samtools_sort/output_unpairedR1
      unpaired_R2:
        source: samtools_sort/output_unpairedR2
      threads: threads_samtools
    scatter:
      - paired
      - unpaired_R1
      - unpaired_R2
    scatterMethod: dotproduct
    out: [output]

  samtools_index:
    run: ../abstract_operations/samtools_index.cwl
    in:
      input:
        source: samtools_merge/output
    scatter: input
    out: [output]

  picard_markduplicates:
    run: ../abstract_operations/subworkflows/picard_markduplicates.cwl
    in:
      input:
        source: samtools_index/output
    scatter: input
    out: [alignments, metrics]

  bam_filtering:
    run: ../abstract_operations/subworkflows/bam_filtering.cwl
    in:
      input:
        source: picard_markduplicates/alignments
      min_mapping_quality: min_mapping_quality
      bits_set: bits_set
      threads: threads_samtools
    scatter: input
    out: [output]

  cnv_manta:
    run: ../abstract_operations/subworkflows/cnv_manta.cwl
    in:
      bams:
        source: bam_filtering/output
      samples: samples
      reference_genome:
        source: bwa_index/output
      bed: bed
      exome: manta_exome
      min_len: manta_min_len
      max_len: manta_max_len
      min_q: manta_min_q
      enable_manta: enable_manta
    when: $(inputs.enable_manta)
    out: [output]

  cnv_gridss:
    run: ../abstract_operations/subworkflows/cnv_gridss.cwl
    in:
      bams:
        source: bam_filtering/output
      samples: samples
      reference_genome:
        source: bwa_index/output
      blacklist: blacklist
      threads: threads_gridss
      min_len: gridss_min_len
      max_len: gridss_max_len
      min_q: gridss_min_q
      enable_gridss: enable_gridss
    when: $(inputs.enable_gridss)
    out: [output]

  cnv_exomedepth:
    run: ../abstract_operations/subworkflows/cnv_exome_depth.cwl
    in:
      bams:
        source: bam_filtering/output
      samples: samples
      reference_genome:
        source: bwa_index/output
      min_len: exomeDepth_min_len
      max_len: exomeDepth_max_len
      min_bf: exomeDepth_min_bf
      enable_exomeDepth: enable_exomeDepth
    when: $(inputs.enable_exomeDepth)
    out: [output]

  cnv_codex:
    run: ../abstract_operations/subworkflows/cnv_codex.cwl
    in:
      bams:
        source: bam_filtering/output
      samples: samples
      bed: bed
      chromosome: chromosome
      min_len: codex_min_len
      max_len: codex_max_len
      min_lratio: codex_min_lratio
      enable_codex: enable_codex
    when: $(inputs.enable_codex)
    out: [output]

  final_filtering:
    run: ../abstract_operations/subworkflows/final_filtering.cwl
    in:
      manta_input:
        source: cnv_manta/output
      gridss_input:
        source: cnv_gridss/output
      exomeDepth_input:
        source: cnv_exomedepth/output
      codex_input:
        source: cnv_codex/output
    out: [output]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0003-4929-1219
    s:email: mailto:laura.rodriguez@bsc.es
    s:name: Laura Rodríguez-Navas
  - class: s:Person
    s:identifier: https://orcid.org/0000-0001-7672-7335
    s:email: mailto:daniel.lopez.lopez@juntadeandalucia.es
    s:name: Daniel López-López

s:codeRepository: https://gitlab.bsc.es/lrodrig1/structuralvariants_poc.git
s:dateCreated: "2020-10-08"
s:license: https://spdx.org/licenses/Apache-2.0
