nextflow.enable.dsl=2

include { COMBINE } from '../modules/combine'
include { HMMSEARCH } from '../modules/hmmsearch'
include { MAKE_TAB_SEP } from '../modules/make_tab_sep'
include { SPLIT_SEQS } from '../modules/split_seqs'

workflow RUN_HMMER {

    take:
    ch_cgc_predicted_proteins
    ch_hmm_database
    ch_hmm_gathering_bit_score
    ch_hmm_omit_alignment
    ch_chunk_size
    ch_name_hmmer
    ch_previous_step_result

    main:
    COMBINE(
        HMMSEARCH.out.output_table.toList(),
        ch_cgc_predicted_proteins,
        ch_name_hmmer
    )

    HMMSEARCH(
        SPLIT_SEQS.out.chunks.flatten().first(),
        ch_hmm_database,
        ch_hmm_gathering_bit_score,
        ch_hmm_omit_alignment
    )

    MAKE_TAB_SEP(
        COMBINE.out.result
    )

    SPLIT_SEQS(
        ch_cgc_predicted_proteins,
        ch_chunk_size
    )

    emit:
    hmm_result = MAKE_TAB_SEP.out.output_with_tabs

}
