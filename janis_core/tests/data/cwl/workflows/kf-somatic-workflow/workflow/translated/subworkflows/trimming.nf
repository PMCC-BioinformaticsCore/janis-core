nextflow.enable.dsl=2

include { COUNT_OVERLAPPED_READS } from '../modules/count_overlapped_reads'
include { TOUCH_EMPTY_FASTA } from '../modules/touch_empty_fasta'
include { TRIMMING } from './trimming'

workflow TRIMMING {

    take:
    ch_reads

    main:
    COUNT_OVERLAPPED_READS(
        ch_reads
    )

    TOUCH_EMPTY_FASTA()

    TRIMMING(
        ch_reads
    )

    emit:
    trimmed_and_reformatted_reads = [TRIMMING.out.trimmed_and_reformatted_reads, TOUCH_EMPTY_FASTA.out.created_file]

}
