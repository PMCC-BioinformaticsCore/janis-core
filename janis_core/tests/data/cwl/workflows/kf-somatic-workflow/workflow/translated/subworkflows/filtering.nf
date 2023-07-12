nextflow.enable.dsl=2

include { COUNT_READS } from '../modules/count_reads'
include { COUNT_READS_AFTER_FILTERING } from '../modules/count_reads_after_filtering'
include { FILTER_CONTIGS_ANTISMASH } from '../modules/filter_contigs_antismash'

workflow FILTERING {

    take:
    ch_contig_min_limit
    ch_fasta

    main:
    COUNT_READS(
        ch_fasta
    )

    COUNT_READS_AFTER_FILTERING(
        FILTER_CONTIGS_ANTISMASH.out.filtered_file
    )

    FILTER_CONTIGS_ANTISMASH(
        ch_fasta,
        ch_contig_min_limit,
        COUNT_READS.out.count
    )

    emit:
    count_after_filtering = COUNT_READS_AFTER_FILTERING.out.count
    filtered_fasta_for_antismash = FILTER_CONTIGS_ANTISMASH.out.filtered_file

}
