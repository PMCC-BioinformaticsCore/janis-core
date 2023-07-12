nextflow.enable.dsl=2

include { CLASSIFICATIONS_TO_OTU_COUNTS } from '../modules/classifications_to_otu_counts'
include { COMPRESS_MAPSEQ } from '../modules/compress_mapseq'
include { COUNT_LINES_MAPSEQ } from '../modules/count_lines_mapseq'
include { COUNTS_TO_HDF5 } from '../modules/counts_to_hdf5'
include { COUNTS_TO_JSON } from '../modules/counts_to_json'
include { MAPSEQ } from '../modules/mapseq'
include { RETURN_OUTPUT_DIR } from '../modules/return_output_dir'
include { VISUALIZE_OTU_COUNTS } from '../modules/visualize_otu_counts'

workflow RUN_ITSONEDB {

    take:
    ch_fasta
    ch_file_for_prefix
    ch_mapseq_ref
    ch_mapseq_taxonomy
    ch_otu_label
    ch_otu_ref
    ch_return_dirname

    main:
    CLASSIFICATIONS_TO_OTU_COUNTS(
        MAPSEQ.out.classifications,
        ch_otu_label,
        ch_otu_ref
    )

    COMPRESS_MAPSEQ(
        MAPSEQ.out.classifications
    )

    COUNT_LINES_MAPSEQ(
        MAPSEQ.out.classifications
    )

    COUNTS_TO_HDF5(
        CLASSIFICATIONS_TO_OTU_COUNTS.out.otu_tsv_notaxid
    )

    COUNTS_TO_JSON(
        CLASSIFICATIONS_TO_OTU_COUNTS.out.otu_tsv_notaxid
    )

    MAPSEQ(
        ch_mapseq_ref,
        ch_file_for_prefix,
        ch_fasta,
        ch_mapseq_taxonomy
    )

    RETURN_OUTPUT_DIR(
        COMPRESS_MAPSEQ.out.compressed_file.toList(),
        ch_return_dirname
    )

    VISUALIZE_OTU_COUNTS(
        CLASSIFICATIONS_TO_OTU_COUNTS.out.otu_txt
    )

    emit:
    number_lines_mapseq = COUNT_LINES_MAPSEQ.out.number
    out_dir = RETURN_OUTPUT_DIR.out.out

}
