nextflow.enable.dsl=2

include { EXTRACT_COORDS } from '../modules/extract_coords'
include { GET_COORDS } from '../modules/get_coords'
include { GET_NCRNAS } from '../modules/get_ncrnas'
include { GZIP_FILES } from '../modules/gzip_files'
include { INDEX_READS } from '../modules/index_reads'
include { RENAME_NCRNAS } from '../modules/rename_ncrnas'

workflow OTHER_NCRNAS {

    take:
    ch_cmsearch_file
    ch_input_sequences
    ch_name_string
    ch_other_ncrna_ribosomal_models

    main:
    EXTRACT_COORDS(
        ch_cmsearch_file,
        ch_name_string
    )

    GET_COORDS(
        EXTRACT_COORDS.out.matched_seqs_with_coords,
        ch_other_ncrna_ribosomal_models
    )

    GET_NCRNAS(
        INDEX_READS.out.sequences_with_index,
        GET_COORDS.out.matches.flatten().first()
    )

    GZIP_FILES(
        RENAME_NCRNAS.out.renamed_file
    )

    INDEX_READS(
        ch_input_sequences
    )

    RENAME_NCRNAS(
        GET_NCRNAS.out.sequences
    )

    emit:
    ncrnas = GZIP_FILES.out.compressed_file

}


workflow OTHER_NCRNAS {

    take:
    ch_cmsearch_file
    ch_input_sequences
    ch_name_string
    ch_other_ncrna_ribosomal_models

    main:
    EXTRACT_COORDS(
        ch_cmsearch_file,
        ch_name_string
    )

    EXTRACT_COORDS(
        ch_cmsearch_file,
        ch_name_string
    )

    GET_COORDS(
        EXTRACT_COORDS.out.matched_seqs_with_coords,
        ch_other_ncrna_ribosomal_models
    )

    GET_COORDS(
        EXTRACT_COORDS.out.matched_seqs_with_coords,
        ch_other_ncrna_ribosomal_models
    )

    GET_NCRNAS(
        INDEX_READS.out.sequences_with_index,
        GET_COORDS.out.matches.flatten().first()
    )

    GET_NCRNAS(
        INDEX_READS.out.sequences_with_index,
        GET_COORDS.out.matches.flatten().first()
    )

    GZIP_FILES(
        RENAME_NCRNAS.out.renamed_file
    )

    GZIP_FILES(
        RENAME_NCRNAS.out.renamed_file
    )

    INDEX_READS(
        ch_input_sequences
    )

    INDEX_READS(
        ch_input_sequences
    )

    RENAME_NCRNAS(
        GET_NCRNAS.out.sequences
    )

    RENAME_NCRNAS(
        GET_NCRNAS.out.sequences
    )

    emit:
    ncrnas = GZIP_FILES.out.compressed_file

}
