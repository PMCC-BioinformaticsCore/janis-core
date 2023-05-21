nextflow.enable.dsl=2

include { COUNT_SUBMITTED_READS } from '../modules/count_submitted_reads'
include { COUNT_SUBMITTED_READS_SINGLE } from '../modules/count_submitted_reads_single'
include { FILTER_PAIRED } from '../modules/filter_paired'
include { OVERLAP_READS } from '../modules/overlap_reads'
include { UNZIP_MERGED_READS } from '../modules/unzip_merged_reads'
include { UNZIP_SINGLE_READS } from '../modules/unzip_single_reads'

workflow OVERLAP_READS {

    take:
    ch_paired_reads_length_filter
    ch_forward_reads
    ch_reverse_reads
    ch_single_reads

    main:
    COUNT_SUBMITTED_READS(
        ch_forward_reads
    )

    COUNT_SUBMITTED_READS_SINGLE(
        UNZIP_SINGLE_READS.out.unzipped_file
    )

    FILTER_PAIRED(
        ch_forward_reads,
        ch_reverse_reads,
        ch_paired_reads_length_filter
    )

    OVERLAP_READS(
        FILTER_PAIRED.out.out_fastq1,
        ch_forward_reads,
        FILTER_PAIRED.out.out_fastq2
    )

    UNZIP_MERGED_READS(
        OVERLAP_READS.out.merged_reads
    )

    UNZIP_SINGLE_READS(
        ch_single_reads
    )

    emit:
    count_forward_submitted_reads = [COUNT_SUBMITTED_READS.out.count, COUNT_SUBMITTED_READS_SINGLE.out.count]
    fastp_report = FILTER_PAIRED.out.json_report
    unzipped_single_reads = [UNZIP_MERGED_READS.out.unzipped_file, UNZIP_SINGLE_READS.out.unzipped_file]

}


workflow OVERLAP_READS {

    take:
    ch_paired_reads_length_filter
    ch_forward_reads
    ch_reverse_reads
    ch_single_reads

    main:
    COUNT_SUBMITTED_READS(
        ch_forward_reads
    )

    COUNT_SUBMITTED_READS(
        ch_forward_reads
    )

    COUNT_SUBMITTED_READS_SINGLE(
        UNZIP_SINGLE_READS.out.unzipped_file
    )

    COUNT_SUBMITTED_READS_SINGLE(
        UNZIP_SINGLE_READS.out.unzipped_file
    )

    FILTER_PAIRED(
        ch_forward_reads,
        ch_reverse_reads,
        ch_paired_reads_length_filter
    )

    FILTER_PAIRED(
        ch_forward_reads,
        ch_reverse_reads,
        ch_paired_reads_length_filter
    )

    OVERLAP_READS(
        FILTER_PAIRED.out.out_fastq1,
        ch_forward_reads,
        FILTER_PAIRED.out.out_fastq2
    )

    OVERLAP_READS(
        FILTER_PAIRED.out.out_fastq1,
        ch_forward_reads,
        FILTER_PAIRED.out.out_fastq2
    )

    UNZIP_MERGED_READS(
        OVERLAP_READS.out.merged_reads
    )

    UNZIP_MERGED_READS(
        OVERLAP_READS.out.merged_reads
    )

    UNZIP_SINGLE_READS(
        ch_single_reads
    )

    UNZIP_SINGLE_READS(
        ch_single_reads
    )

    emit:
    count_forward_submitted_reads = [COUNT_SUBMITTED_READS.out.count, COUNT_SUBMITTED_READS_SINGLE.out.count]
    fastp_report = FILTER_PAIRED.out.json_report
    unzipped_single_reads = [UNZIP_MERGED_READS.out.unzipped_file, UNZIP_SINGLE_READS.out.unzipped_file]

}
