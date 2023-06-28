nextflow.enable.dsl=2

include { CLASSIFY_LSUS } from './classify_LSUs'
include { CLASSIFY_SSUS } from './classify_SSUs'
include { COUNT_LSU_FASTA } from '../modules/count_lsu_fasta'
include { COUNT_SSU_FASTA } from '../modules/count_ssu_fasta'
include { EXTRACT_COORDS } from '../modules/extract_coords'
include { EXTRACT_SEQUENCES } from '../modules/extract_sequences'
include { EXTRACT_SUBUNITS } from '../modules/extract_subunits'
include { EXTRACT_SUBUNITS_COORDS } from '../modules/extract_subunits_coords'
include { FIND_RIBOSOMAL_NC_RNAS } from './find_ribosomal_ncRNAs'
include { GZIP_FILES } from '../modules/gzip_files'
include { INDEX_READS } from '../modules/index_reads'

workflow RNA_PREDICTION {

    take:
    ch_input_sequences
    ch_ncrna_ribosomal_model_clans
    ch_ncrna_ribosomal_models
    ch_pattern_5.8s
    ch_pattern_5_s
    ch_pattern_lsu
    ch_pattern_ssu
    ch_silva_lsu_database
    ch_silva_lsu_otus
    ch_silva_lsu_taxonomy
    ch_silva_ssu_database
    ch_silva_ssu_otus
    ch_silva_ssu_taxonomy
    ch_type

    main:
    CLASSIFY_LSUS(
        EXTRACT_SUBUNITS.out.LSU_seqs,
        ch_input_sequences,
        ch_silva_lsu_database,
        ch_silva_lsu_taxonomy,
        ch_pattern_lsu,
        ch_silva_lsu_otus,
        params.after_qc.rna_prediction.classify_lsus_return_dirname
    )

    CLASSIFY_SSUS(
        EXTRACT_SUBUNITS.out.SSU_seqs,
        ch_input_sequences,
        ch_silva_ssu_database,
        ch_silva_ssu_taxonomy,
        ch_pattern_ssu,
        ch_silva_ssu_otus,
        params.after_qc.rna_prediction.classify_ssus_return_dirname
    )

    COUNT_LSU_FASTA(
        EXTRACT_SUBUNITS.out.LSU_seqs
    )

    COUNT_SSU_FASTA(
        EXTRACT_SUBUNITS.out.SSU_seqs
    )

    EXTRACT_COORDS(
        FIND_RIBOSOMAL_NC_RNAS.out.deoverlapped_matches
    )

    EXTRACT_SEQUENCES(
        INDEX_READS.out.sequences_with_index,
        EXTRACT_COORDS.out.matched_seqs_with_coords
    )

    EXTRACT_SUBUNITS(
        EXTRACT_SEQUENCES.out.sequences,
        ch_pattern_lsu,
        ch_pattern_ssu,
        ch_pattern_5.8s,
        ch_pattern_5_s
    )

    EXTRACT_SUBUNITS_COORDS(
        EXTRACT_COORDS.out.matched_seqs_with_coords,
        ch_pattern_lsu,
        ch_pattern_ssu
    )

    FIND_RIBOSOMAL_NC_RNAS(
        ch_ncrna_ribosomal_model_clans,
        ch_ncrna_ribosomal_models,
        ch_input_sequences,
        ch_type
    )

    GZIP_FILES(
        EXTRACT_SUBUNITS.out.fastas.flatten().first()
    )

    INDEX_READS(
        ch_input_sequences
    )

    emit:
    LSU_SSU_count = EXTRACT_SUBUNITS_COORDS.out.counts
    LSU_coords = EXTRACT_SUBUNITS_COORDS.out.LSU_seqs
    LSU_fasta = EXTRACT_SUBUNITS.out.LSU_seqs
    LSU_folder = CLASSIFY_LSUS.out.out_dir
    SSU_coords = EXTRACT_SUBUNITS_COORDS.out.SSU_seqs
    SSU_fasta = EXTRACT_SUBUNITS.out.SSU_seqs
    SSU_folder = CLASSIFY_SSUS.out.out_dir
    cmsearch_result = FIND_RIBOSOMAL_NC_RNAS.out.concatenate_matches
    compressed_rnas = GZIP_FILES.out.compressed_file
    ncRNA = FIND_RIBOSOMAL_NC_RNAS.out.deoverlapped_matches
    number_LSU_mapseq = CLASSIFY_LSUS.out.number_lines_mapseq
    number_SSU_mapseq = CLASSIFY_SSUS.out.number_lines_mapseq

}


workflow RNA_PREDICTION {

    take:
    ch_input_sequences
    ch_ncrna_ribosomal_model_clans
    ch_ncrna_ribosomal_models
    ch_pattern_5.8s
    ch_pattern_5_s
    ch_pattern_lsu
    ch_pattern_ssu
    ch_silva_lsu_database
    ch_silva_lsu_otus
    ch_silva_lsu_taxonomy
    ch_silva_ssu_database
    ch_silva_ssu_otus
    ch_silva_ssu_taxonomy
    ch_type

    main:
    CLASSIFY_LSUS(
        EXTRACT_SUBUNITS.out.LSU_seqs,
        ch_input_sequences,
        ch_silva_lsu_database,
        ch_silva_lsu_taxonomy,
        ch_pattern_lsu,
        ch_silva_lsu_otus,
        params.after_qc.rna_prediction.classify_lsus_return_dirname
    )

    CLASSIFY_LSUS(
        EXTRACT_SUBUNITS.out.LSU_seqs,
        ch_input_sequences,
        ch_silva_lsu_database,
        ch_silva_lsu_taxonomy,
        ch_pattern_lsu,
        ch_silva_lsu_otus,
        params.after_qc.rna_prediction.classify_lsus_return_dirname
    )

    CLASSIFY_SSUS(
        EXTRACT_SUBUNITS.out.SSU_seqs,
        ch_input_sequences,
        ch_silva_ssu_database,
        ch_silva_ssu_taxonomy,
        ch_pattern_ssu,
        ch_silva_ssu_otus,
        params.after_qc.rna_prediction.classify_ssus_return_dirname
    )

    CLASSIFY_SSUS(
        EXTRACT_SUBUNITS.out.SSU_seqs,
        ch_input_sequences,
        ch_silva_ssu_database,
        ch_silva_ssu_taxonomy,
        ch_pattern_ssu,
        ch_silva_ssu_otus,
        params.after_qc.rna_prediction.classify_ssus_return_dirname
    )

    COUNT_LSU_FASTA(
        EXTRACT_SUBUNITS.out.LSU_seqs
    )

    COUNT_LSU_FASTA(
        EXTRACT_SUBUNITS.out.LSU_seqs
    )

    COUNT_SSU_FASTA(
        EXTRACT_SUBUNITS.out.SSU_seqs
    )

    COUNT_SSU_FASTA(
        EXTRACT_SUBUNITS.out.SSU_seqs
    )

    EXTRACT_COORDS(
        FIND_RIBOSOMAL_NC_RNAS.out.deoverlapped_matches
    )

    EXTRACT_COORDS(
        FIND_RIBOSOMAL_NC_RNAS.out.deoverlapped_matches
    )

    EXTRACT_SEQUENCES(
        INDEX_READS.out.sequences_with_index,
        EXTRACT_COORDS.out.matched_seqs_with_coords
    )

    EXTRACT_SEQUENCES(
        INDEX_READS.out.sequences_with_index,
        EXTRACT_COORDS.out.matched_seqs_with_coords
    )

    EXTRACT_SUBUNITS(
        EXTRACT_SEQUENCES.out.sequences,
        ch_pattern_lsu,
        ch_pattern_ssu,
        ch_pattern_5.8s,
        ch_pattern_5_s
    )

    EXTRACT_SUBUNITS(
        EXTRACT_SEQUENCES.out.sequences,
        ch_pattern_lsu,
        ch_pattern_ssu,
        ch_pattern_5.8s,
        ch_pattern_5_s
    )

    EXTRACT_SUBUNITS_COORDS(
        EXTRACT_COORDS.out.matched_seqs_with_coords,
        ch_pattern_lsu,
        ch_pattern_ssu
    )

    EXTRACT_SUBUNITS_COORDS(
        EXTRACT_COORDS.out.matched_seqs_with_coords,
        ch_pattern_lsu,
        ch_pattern_ssu
    )

    FIND_RIBOSOMAL_NC_RNAS(
        ch_ncrna_ribosomal_model_clans,
        ch_ncrna_ribosomal_models,
        ch_input_sequences,
        ch_type
    )

    FIND_RIBOSOMAL_NC_RNAS(
        ch_ncrna_ribosomal_model_clans,
        ch_ncrna_ribosomal_models,
        ch_input_sequences,
        ch_type
    )

    GZIP_FILES(
        EXTRACT_SUBUNITS.out.fastas.flatten().first()
    )

    GZIP_FILES(
        EXTRACT_SUBUNITS.out.fastas.flatten().first()
    )

    INDEX_READS(
        ch_input_sequences
    )

    INDEX_READS(
        ch_input_sequences
    )

    emit:
    LSU_SSU_count = EXTRACT_SUBUNITS_COORDS.out.counts
    LSU_coords = EXTRACT_SUBUNITS_COORDS.out.LSU_seqs
    LSU_fasta = EXTRACT_SUBUNITS.out.LSU_seqs
    LSU_folder = CLASSIFY_LSUS.out.out_dir
    SSU_coords = EXTRACT_SUBUNITS_COORDS.out.SSU_seqs
    SSU_fasta = EXTRACT_SUBUNITS.out.SSU_seqs
    SSU_folder = CLASSIFY_SSUS.out.out_dir
    cmsearch_result = FIND_RIBOSOMAL_NC_RNAS.out.concatenate_matches
    compressed_rnas = GZIP_FILES.out.compressed_file
    ncRNA = FIND_RIBOSOMAL_NC_RNAS.out.deoverlapped_matches
    number_LSU_mapseq = CLASSIFY_LSUS.out.number_lines_mapseq
    number_SSU_mapseq = CLASSIFY_SSUS.out.number_lines_mapseq

}
