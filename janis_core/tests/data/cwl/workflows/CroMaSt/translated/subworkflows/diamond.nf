nextflow.enable.dsl=2

include { DIAMOND_RUN } from '../modules/diamond_run'
include { POST_PROCESSING_UNIREF90 } from './post_processing_uniref90'

workflow DIAMOND {

    take:
    ch_uniref90_db_txt
    ch_database_file
    ch_filename
    ch_max_target_seqs
    ch_output_format
    ch_query_input_file
    ch_strand
    ch_threads

    main:
    DIAMOND_RUN(
        ch_query_input_file,
        ch_database_file,
        ch_max_target_seqs,
        ch_output_format,
        ch_strand,
        ch_threads
    )

    POST_PROCESSING_UNIREF90(
        ch_filename,
        ch_uniref90_db_txt,
        DIAMOND_RUN.out.matches
    )

    emit:
    diamond_output = DIAMOND_RUN.out.matches
    post_processing_output = POST_PROCESSING_UNIREF90.out.join_out

}
