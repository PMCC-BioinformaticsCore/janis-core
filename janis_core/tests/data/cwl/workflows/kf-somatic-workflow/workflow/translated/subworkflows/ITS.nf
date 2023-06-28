nextflow.enable.dsl=2

include { CAT } from '../modules/cat'
include { COUNT_ITS_SEQS } from '../modules/count_ITS_seqs'
include { COUNT_MASKED_FASTA } from '../modules/count_masked_fasta'
include { GZIP_MASKED_ITS } from '../modules/gzip_masked_ITS'
include { MASK_FOR_ITS } from '../modules/mask_for_ITS'
include { REFORMAT_COORDS } from '../modules/reformat_coords'
include { RUN_ITSONEDB } from './run_itsonedb'
include { RUN_UNITE } from './run_unite'

workflow ITS {

    take:
    ch_lsu_coordinates
    ch_ssu_coordinates
    ch_itsone_database
    ch_itsone_otus
    ch_itsone_taxonomy
    ch_otu_itsone_label
    ch_otu_unite_label
    ch_query_sequences
    ch_unite_database
    ch_unite_otus
    ch_unite_taxonomy

    main:
    CAT(
        ch_ssu_coordinates.toList()
    )

    COUNT_ITS_SEQS(
        MASK_FOR_ITS.out.masked_sequences
    )

    COUNT_MASKED_FASTA(
        MASK_FOR_ITS.out.masked_sequences
    )

    GZIP_MASKED_ITS(
        MASK_FOR_ITS.out.masked_sequences
    )

    MASK_FOR_ITS(
        REFORMAT_COORDS.out.maskfile,
        ch_query_sequences
    )

    REFORMAT_COORDS(
        CAT.out.result
    )

    RUN_ITSONEDB(
        MASK_FOR_ITS.out.masked_sequences,
        ch_query_sequences,
        ch_itsone_database,
        ch_itsone_taxonomy,
        ch_otu_itsone_label,
        ch_itsone_otus,
        params.after_qc.its.run_itsonedb_return_dirname
    )

    RUN_UNITE(
        MASK_FOR_ITS.out.masked_sequences,
        ch_query_sequences,
        ch_unite_database,
        ch_unite_taxonomy,
        ch_otu_unite_label,
        ch_unite_otus,
        params.after_qc.its.run_unite_return_dirname
    )

    emit:
    itsonedb_folder = RUN_ITSONEDB.out.out_dir
    masking_file = GZIP_MASKED_ITS.out.compressed_file
    number_ITS_seqs = COUNT_ITS_SEQS.out.count
    unite_folder = RUN_UNITE.out.out_dir

}
