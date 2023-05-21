nextflow.enable.dsl=2

include { QC_FLAG } from '../modules/QC_FLAG'
include { COUNT_PROCESSED_READS } from '../modules/count_processed_reads'
include { HASHSUM_PAIRED } from '../modules/hashsum_paired'
include { HASHSUM_SINGLE } from '../modules/hashsum_single'
include { OVERLAP_READS } from './overlap_reads'
include { QC_STATS } from '../modules/qc_stats'
include { RUN_QUALITY_CONTROL_FILTERING } from '../modules/run_quality_control_filtering'
include { TRIMMING } from './trimming'
include { CLEAN_FASTA_HEADERS } from '../modules/clean_fasta_headers'
include { CONVERT_TRIMMED_READS_TO_FASTA } from '../modules/convert_trimmed_reads_to_fasta'
include { LENGTH_FILTER } from '../modules/length_filter'
include { TRIM_QUALITY_CONTROL } from '../modules/trim_quality_control'

workflow BEFORE_QC {

    take:
    ch_qc_min_length
    ch_stats_file_name
    ch_forward_reads
    ch_reverse_reads
    ch_single_reads

    main:
    QC_FLAG(
        COUNT_PROCESSED_READS.out.count
    )

    COUNT_PROCESSED_READS(
        RUN_QUALITY_CONTROL_FILTERING.out.filtered_file
    )

    HASHSUM_PAIRED(
        ch_forward_reads
    )

    HASHSUM_SINGLE(
        ch_single_reads
    )

    OVERLAP_READS(
        params.before_qc.overlap_reads_paired_reads_length_filter,
        ch_forward_reads,
        ch_reverse_reads,
        ch_single_reads
    )

    QC_STATS(
        RUN_QUALITY_CONTROL_FILTERING.out.filtered_file,
        COUNT_PROCESSED_READS.out.count
    )

    RUN_QUALITY_CONTROL_FILTERING(
        TRIMMING.out.trimmed_and_reformatted_reads,
        ch_qc_min_length,
        OVERLAP_READS.out.count_forward_submitted_reads
    )

    TRIMMING(
        OVERLAP_READS.out.unzipped_single_reads
    )

    emit:
    fastp_filtering_json = OVERLAP_READS.out.fastp_report
    filtered_fasta = RUN_QUALITY_CONTROL_FILTERING.out.filtered_file
    input_files_hashsum_paired = HASHSUM_PAIRED.out.hashsum
    input_files_hashsum_single = HASHSUM_SINGLE.out.hashsum
    qc_statistics = QC_STATS.out.output_dir
    qc_status = QC_FLAG.out.qc_flag
    qc_summary = RUN_QUALITY_CONTROL_FILTERING.out.stats_summary_file

}


workflow BEFORE_QC {

    take:
    ch_qc_min_length
    ch_forward_reads
    ch_reverse_reads
    ch_single_reads

    main:
    QC_FLAG(
        COUNT_PROCESSED_READS.out.count
    )

    QC_FLAG(
        COUNT_PROCESSED_READS.out.count
    )

    CLEAN_FASTA_HEADERS(
        TRIM_QUALITY_CONTROL.out.reads1_trimmed
    )

    CONVERT_TRIMMED_READS_TO_FASTA(
        CLEAN_FASTA_HEADERS.out.sequences_with_cleaned_headers
    )

    COUNT_PROCESSED_READS(
        LENGTH_FILTER.out.filtered_file
    )

    COUNT_PROCESSED_READS(
        LENGTH_FILTER.out.filtered_file
    )

    HASHSUM_PAIRED(
        ch_forward_reads
    )

    HASHSUM_PAIRED(
        ch_forward_reads
    )

    HASHSUM_SINGLE(
        ch_single_reads
    )

    HASHSUM_SINGLE(
        ch_single_reads
    )

    LENGTH_FILTER(
        CONVERT_TRIMMED_READS_TO_FASTA.out.fasta,
        ch_qc_min_length,
        OVERLAP_READS.out.count_forward_submitted_reads
    )

    OVERLAP_READS(
        params.before_qc.overlap_reads_paired_reads_length_filter,
        ch_forward_reads,
        ch_reverse_reads,
        ch_single_reads
    )

    OVERLAP_READS(
        params.before_qc.overlap_reads_paired_reads_length_filter,
        ch_forward_reads,
        ch_reverse_reads,
        ch_single_reads
    )

    QC_STATS(
        LENGTH_FILTER.out.filtered_file,
        COUNT_PROCESSED_READS.out.count
    )

    QC_STATS(
        LENGTH_FILTER.out.filtered_file,
        COUNT_PROCESSED_READS.out.count
    )

    TRIM_QUALITY_CONTROL(
        OVERLAP_READS.out.unzipped_single_reads
    )

    emit:
    fastp_filtering_json = OVERLAP_READS.out.fastp_report
    filtered_fasta = LENGTH_FILTER.out.filtered_file
    input_files_hashsum_paired = HASHSUM_PAIRED.out.hashsum
    input_files_hashsum_single = HASHSUM_SINGLE.out.hashsum
    motus_input = CLEAN_FASTA_HEADERS.out.sequences_with_cleaned_headers
    qc_statistics = QC_STATS.out.output_dir
    qc_status = QC_FLAG.out.qc_flag
    qc_summary = LENGTH_FILTER.out.stats_summary_file

}
