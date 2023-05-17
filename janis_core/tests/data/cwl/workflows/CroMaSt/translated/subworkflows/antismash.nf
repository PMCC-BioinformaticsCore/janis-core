nextflow.enable.dsl=2

include { CHUNKING } from './chunking'
include { FILTERING } from './filtering'
include { NO_ANTISMASH_SUBWF } from './no_antismash_subwf'

workflow ANTISMASH {

    take:
    ch_clusters_glossary
    ch_final_folder_name
    ch_input_filtered_fasta

    main:
    CHUNKING(
        ch_clusters_glossary,
        FILTERING.out.filtered_fasta_for_antismash,
        ch_final_folder_name,
        params.after_qc.antismash.chunking_split_size
    )

    FILTERING(
        params.after_qc.antismash.contig_min_limit,
        ch_input_filtered_fasta
    )

    NO_ANTISMASH_SUBWF(
        FILTERING.out.filtered_fasta_for_antismash,
        ch_final_folder_name
    )

    emit:
    antismash_clusters = CHUNKING.out.antismash_clusters
    antismash_folder = [NO_ANTISMASH_SUBWF.out.antismash_result_folder, CHUNKING.out.antismash_folder_chunking]

}
