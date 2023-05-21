nextflow.enable.dsl=2

include { JOIN } from '../modules/join'
include { SORTING } from '../modules/sorting'

workflow POST_PROCESSING_UNIREF90 {

    take:
    ch_filename
    ch_input_db
    ch_input_diamond

    main:
    JOIN(
        SORTING.out.output_sorted,
        ch_filename,
        ch_input_db
    )

    SORTING(
        ch_input_diamond
    )

    emit:
    join_out = JOIN.out.output_join

}
