nextflow.enable.dsl=2

include { ITS } from './ITS'
include { GZIP_LSU } from '../modules/gzip_LSU'
include { GZIP_SSU } from '../modules/gzip_SSU'
include { GZIP_FILES } from '../modules/gzip_files'
include { NO_TAX_FILE_FLAG } from '../modules/no_tax_file_flag'
include { RETURN_ITS_DIR } from '../modules/return_its_dir'
include { RETURN_SEQ_DIR } from '../modules/return_seq_dir'
include { RNA_PREDICTION } from './rna_prediction'
include { SUPPRESS_TAX } from '../modules/suppress_tax'
include { ANTISMASH } from './antismash'
include { CGC } from './cgc'
include { CHUNKING_FINAL } from './chunking_final'
include { COMPRESSION } from '../modules/compression'
include { FASTA_INDEX } from '../modules/fasta_index'
include { FUNCTIONAL_ANNOTATION_AND_POST_PROCESSING } from './functional_annotation_and_post_processing'
include { MOVE_TO_SEQ_CAT_FOLDER } from '../modules/move_to_seq_cat_folder'
include { OTHER_NCRNAS } from './other_ncrnas'
include { RETURN_TAX_DIR } from '../modules/return_tax_dir'
include { MOTUS_TAXONOMY } from './motus_taxonomy'

workflow AFTER_QC {

    take:
    ch_5.8s_pattern
    ch_5s_pattern
    ch_filtered_fasta
    ch_itsonedb
    ch_itsonedb_label
    ch_itsonedb_otu_file
    ch_itsonedb_tax
    ch_lsu_db
    ch_lsu_label
    ch_lsu_otus
    ch_lsu_tax
    ch_rfam_model_clans
    ch_rfam_models
    ch_ssu_db
    ch_ssu_label
    ch_ssu_otus
    ch_ssu_tax
    ch_unite_db
    ch_unite_label
    ch_unite_otu_file
    ch_unite_tax

    main:
    ITS(
        RNA_PREDICTION.out.LSU_coords,
        RNA_PREDICTION.out.SSU_coords,
        ch_itsonedb,
        ch_itsonedb_otu_file,
        ch_itsonedb_tax,
        ch_itsonedb_label,
        ch_unite_label,
        ch_filtered_fasta,
        ch_unite_db,
        ch_unite_otu_file,
        ch_unite_tax
    )

    GZIP_LSU(
        RNA_PREDICTION.out.LSU_fasta
    )

    GZIP_SSU(
        RNA_PREDICTION.out.SSU_fasta
    )

    GZIP_FILES(
        ch_filtered_fasta
    )

    NO_TAX_FILE_FLAG()

    RETURN_ITS_DIR(
        ITS.out.unite_folder.toList()
    )

    RETURN_SEQ_DIR(
        RNA_PREDICTION.out.compressed_rnas
    )

    RNA_PREDICTION(
        ch_filtered_fasta,
        ch_rfam_model_clans,
        ch_rfam_models,
        ch_5.8s_pattern,
        ch_5s_pattern,
        ch_lsu_label,
        ch_ssu_label,
        ch_lsu_db,
        ch_lsu_otus,
        ch_lsu_tax,
        ch_ssu_db,
        ch_ssu_otus,
        ch_ssu_tax,
        params.after_qc.rna_prediction_type
    )

    SUPPRESS_TAX(
        ITS.out.masking_file,
        GZIP_LSU.out.compressed_file,
        GZIP_SSU.out.compressed_file,
        RETURN_ITS_DIR.out.out,
        RNA_PREDICTION.out.LSU_folder,
        RNA_PREDICTION.out.SSU_folder
    )

    emit:
    ITS_length = SUPPRESS_TAX.out.its_length
    gz_files = GZIP_FILES.out.compressed_file
    optional_tax_file_flag = NO_TAX_FILE_FLAG.out.created_file
    rna_count = RNA_PREDICTION.out.LSU_SSU_count
    sequence_categorisation_folder = RETURN_SEQ_DIR.out.out
    suppressed_upload = SUPPRESS_TAX.out.out_suppress
    taxonomy_summary_folder = SUPPRESS_TAX.out.out_tax

}


workflow AFTER_QC {

    take:
    ch_5.8s_pattern
    ch_5s_pattern
    ch_cgc_config
    ch_cgc_postfixes
    ch_egg_nog_data_dir
    ch_egg_nog_db
    ch_egg_nog_diamond_db
    ch_hmm_gathering_bit_score
    ch_hmm_name_database
    ch_hmm_omit_alignment
    ch_inter_pro_scan_applications
    ch_inter_pro_scan_databases
    ch_inter_pro_scan_output_format
    ch_uniref90_db_txt
    ch_cgc_chunk_size
    ch_clusters_glossary
    ch_diamond_database_file
    ch_diamond_header
    ch_diamond_max_target_seqs
    ch_filtered_fasta
    ch_func_ann_names_hmmer
    ch_func_ann_names_ips
    ch_go_config
    ch_gp_flatfiles_path
    ch_graphs
    ch_hmmsearch_header
    ch_ips_header
    ch_ko_file
    ch_lsu_db
    ch_lsu_label
    ch_lsu_otus
    ch_lsu_tax
    ch_other_ncrna_models
    ch_pathways_classes
    ch_pathways_names
    ch_protein_chunk_size_ips
    ch_protein_chunk_size_eggnog
    ch_protein_chunk_size_hmm
    ch_rfam_model_clans
    ch_rfam_models
    ch_ssu_db
    ch_ssu_label
    ch_ssu_otus
    ch_ssu_tax

    main:
    ANTISMASH(
        ch_clusters_glossary,
        params.after_qc.antismash_final_folder_name,
        ch_filtered_fasta
    )

    CGC(
        ch_cgc_chunk_size,
        ch_filtered_fasta,
        RNA_PREDICTION.out.ncRNA,
        ch_cgc_postfixes
    )

    CHUNKING_FINAL(
        CGC.out.results.flatten().first(),
        ch_filtered_fasta,
        CGC.out.results.flatten().first(),
        RNA_PREDICTION.out.LSU_fasta,
        RNA_PREDICTION.out.SSU_fasta
    )

    COMPRESSION(
        RNA_PREDICTION.out.ncRNA
    )

    FASTA_INDEX(
        ch_filtered_fasta
    )

    FUNCTIONAL_ANNOTATION_AND_POST_PROCESSING(
        ch_egg_nog_data_dir,
        ch_egg_nog_db,
        ch_egg_nog_diamond_db,
        ch_hmm_gathering_bit_score,
        ch_hmm_name_database,
        ch_hmm_omit_alignment,
        ch_inter_pro_scan_applications,
        ch_inter_pro_scan_databases,
        ch_inter_pro_scan_output_format,
        ch_uniref90_db_txt,
        ANTISMASH.out.antismash_clusters,
        CGC.out.results.flatten().first(),
        ch_diamond_database_file,
        ch_diamond_header,
        ch_diamond_max_target_seqs,
        ch_filtered_fasta,
        ch_func_ann_names_hmmer,
        ch_func_ann_names_ips,
        ch_go_config,
        ch_gp_flatfiles_path,
        ch_graphs,
        ch_hmmsearch_header,
        ch_ips_header,
        ch_ko_file,
        ch_pathways_classes,
        ch_pathways_names,
        ch_protein_chunk_size_ips,
        ch_protein_chunk_size_eggnog,
        ch_protein_chunk_size_hmm,
        RNA_PREDICTION.out.ncRNA
    )

    MOVE_TO_SEQ_CAT_FOLDER(
        CHUNKING_FINAL.out.SC_fasta_chunks
    )

    NO_TAX_FILE_FLAG()

    NO_TAX_FILE_FLAG()

    OTHER_NCRNAS(
        RNA_PREDICTION.out.ncRNA,
        ch_filtered_fasta,
        params.after_qc.other_ncrnas_name_string,
        ch_other_ncrna_models
    )

    RETURN_TAX_DIR(
        RNA_PREDICTION.out.SSU_folder.toList()
    )

    RNA_PREDICTION(
        ch_filtered_fasta,
        ch_rfam_model_clans,
        ch_rfam_models,
        ch_5.8s_pattern,
        ch_5s_pattern,
        ch_lsu_label,
        ch_ssu_label,
        ch_lsu_db,
        ch_lsu_otus,
        ch_lsu_tax,
        ch_ssu_db,
        ch_ssu_otus,
        ch_ssu_tax,
        params.after_qc.rna_prediction_type
    )

    RNA_PREDICTION(
        ch_filtered_fasta,
        ch_rfam_model_clans,
        ch_rfam_models,
        ch_5.8s_pattern,
        ch_5s_pattern,
        ch_lsu_label,
        ch_ssu_label,
        ch_lsu_db,
        ch_lsu_otus,
        ch_lsu_tax,
        ch_ssu_db,
        ch_ssu_otus,
        ch_ssu_tax,
        params.after_qc.rna_prediction_type
    )

    emit:
    bgzip_fasta_file = FASTA_INDEX.out.fasta_bgz
    bgzip_index = FASTA_INDEX.out.bgz_index
    chunking_nucleotides = CHUNKING_FINAL.out.nucleotide_fasta_chunks
    chunking_proteins = CHUNKING_FINAL.out.protein_fasta_chunks
    compressed_files = COMPRESSION.out.compressed_file
    count_CDS = CGC.out.count_faa
    functional_annotation_folder = FUNCTIONAL_ANNOTATION_AND_POST_PROCESSING.out.functional_annotation_folder
    index_fasta_file = FASTA_INDEX.out.fasta_index
    optional_tax_file_flag = NO_TAX_FILE_FLAG.out.created_file
    pathways_systems_folder = FUNCTIONAL_ANNOTATION_AND_POST_PROCESSING.out.pathways_systems_folder
    pathways_systems_folder_antismash = ANTISMASH.out.antismash_folder
    pathways_systems_folder_antismash_summary = FUNCTIONAL_ANNOTATION_AND_POST_PROCESSING.out.pathways_systems_folder_antismash_summary
    rna_count = RNA_PREDICTION.out.LSU_SSU_count
    sequence_categorisation_folder = MOVE_TO_SEQ_CAT_FOLDER.out.out
    stats = FUNCTIONAL_ANNOTATION_AND_POST_PROCESSING.out.stats
    taxonomy_summary_folder = RETURN_TAX_DIR.out.out

}


workflow AFTER_QC {

    take:
    ch_5.8s_pattern
    ch_5s_pattern
    ch_cgc_config
    ch_cgc_postfixes
    ch_hmm_gathering_bit_score
    ch_hmm_name_database
    ch_hmm_omit_alignment
    ch_inter_pro_scan_applications
    ch_inter_pro_scan_databases
    ch_inter_pro_scan_output_format
    ch_cgc_chunk_size
    ch_filtered_fasta
    ch_func_ann_names_hmmer
    ch_func_ann_names_ips
    ch_go_config
    ch_hmmsearch_header
    ch_ips_header
    ch_ko_file
    ch_lsu_db
    ch_lsu_label
    ch_lsu_otus
    ch_lsu_tax
    ch_motus_input
    ch_other_ncrna_models
    ch_protein_chunk_size_ips
    ch_protein_chunk_size_hmm
    ch_rfam_model_clans
    ch_rfam_models
    ch_ssu_db
    ch_ssu_label
    ch_ssu_otus
    ch_ssu_tax
    ch_egg_nog_data_dir
    ch_egg_nog_db
    ch_egg_nog_diamond_db

    main:
    CGC(
        ch_cgc_chunk_size,
        ch_filtered_fasta,
        RNA_PREDICTION.out.ncRNA,
        ch_cgc_postfixes
    )

    CGC(
        ch_cgc_chunk_size,
        ch_filtered_fasta,
        RNA_PREDICTION.out.ncRNA,
        ch_cgc_postfixes
    )

    CHUNKING_FINAL(
        CGC.out.results.flatten().first(),
        ch_filtered_fasta,
        CGC.out.results.flatten().first(),
        RNA_PREDICTION.out.LSU_fasta,
        RNA_PREDICTION.out.SSU_fasta
    )

    CHUNKING_FINAL(
        CGC.out.results.flatten().first(),
        ch_filtered_fasta,
        CGC.out.results.flatten().first(),
        RNA_PREDICTION.out.LSU_fasta,
        RNA_PREDICTION.out.SSU_fasta
    )

    COMPRESSION(
        RNA_PREDICTION.out.ncRNA
    )

    COMPRESSION(
        RNA_PREDICTION.out.ncRNA
    )

    FUNCTIONAL_ANNOTATION_AND_POST_PROCESSING(
        ch_hmm_name_database,
        ch_hmm_gathering_bit_score,
        ch_hmm_omit_alignment,
        ch_inter_pro_scan_applications,
        ch_inter_pro_scan_databases,
        ch_inter_pro_scan_output_format,
        CGC.out.results.flatten().first(),
        ch_filtered_fasta,
        ch_func_ann_names_hmmer,
        ch_func_ann_names_ips,
        ch_go_config,
        ch_hmmsearch_header,
        ch_ips_header,
        ch_ko_file,
        ch_protein_chunk_size_ips,
        ch_protein_chunk_size_hmm,
        RNA_PREDICTION.out.ncRNA
    )

    FUNCTIONAL_ANNOTATION_AND_POST_PROCESSING(
        ch_hmm_name_database,
        ch_hmm_gathering_bit_score,
        ch_hmm_omit_alignment,
        ch_inter_pro_scan_applications,
        ch_inter_pro_scan_databases,
        ch_inter_pro_scan_output_format,
        CGC.out.results.flatten().first(),
        ch_filtered_fasta,
        ch_func_ann_names_hmmer,
        ch_func_ann_names_ips,
        ch_go_config,
        ch_hmmsearch_header,
        ch_ips_header,
        ch_ko_file,
        ch_protein_chunk_size_ips,
        ch_protein_chunk_size_hmm,
        RNA_PREDICTION.out.ncRNA
    )

    MOTUS_TAXONOMY(
        ch_motus_input
    )

    MOVE_TO_SEQ_CAT_FOLDER(
        CHUNKING_FINAL.out.SC_fasta_chunks
    )

    MOVE_TO_SEQ_CAT_FOLDER(
        CHUNKING_FINAL.out.SC_fasta_chunks
    )

    NO_TAX_FILE_FLAG()

    NO_TAX_FILE_FLAG()

    NO_TAX_FILE_FLAG()

    OTHER_NCRNAS(
        RNA_PREDICTION.out.ncRNA,
        ch_filtered_fasta,
        params.after_qc.other_ncrnas_name_string,
        ch_other_ncrna_models
    )

    OTHER_NCRNAS(
        RNA_PREDICTION.out.ncRNA,
        ch_filtered_fasta,
        params.after_qc.other_ncrnas_name_string,
        ch_other_ncrna_models
    )

    RETURN_TAX_DIR(
        RNA_PREDICTION.out.SSU_folder.toList()
    )

    RETURN_TAX_DIR(
        RNA_PREDICTION.out.SSU_folder.toList()
    )

    RNA_PREDICTION(
        ch_filtered_fasta,
        ch_rfam_model_clans,
        ch_rfam_models,
        ch_5.8s_pattern,
        ch_5s_pattern,
        ch_lsu_label,
        ch_ssu_label,
        ch_lsu_db,
        ch_lsu_otus,
        ch_lsu_tax,
        ch_ssu_db,
        ch_ssu_otus,
        ch_ssu_tax,
        params.after_qc.rna_prediction_type
    )

    RNA_PREDICTION(
        ch_filtered_fasta,
        ch_rfam_model_clans,
        ch_rfam_models,
        ch_5.8s_pattern,
        ch_5s_pattern,
        ch_lsu_label,
        ch_ssu_label,
        ch_lsu_db,
        ch_lsu_otus,
        ch_lsu_tax,
        ch_ssu_db,
        ch_ssu_otus,
        ch_ssu_tax,
        params.after_qc.rna_prediction_type
    )

    RNA_PREDICTION(
        ch_filtered_fasta,
        ch_rfam_model_clans,
        ch_rfam_models,
        ch_5.8s_pattern,
        ch_5s_pattern,
        ch_lsu_label,
        ch_ssu_label,
        ch_lsu_db,
        ch_lsu_otus,
        ch_lsu_tax,
        ch_ssu_db,
        ch_ssu_otus,
        ch_ssu_tax,
        params.after_qc.rna_prediction_type
    )

    emit:
    chunking_nucleotides = CHUNKING_FINAL.out.nucleotide_fasta_chunks
    chunking_proteins = CHUNKING_FINAL.out.protein_fasta_chunks
    compressed_files = COMPRESSION.out.compressed_file
    count_CDS = CGC.out.count_faa
    functional_annotation_folder = FUNCTIONAL_ANNOTATION_AND_POST_PROCESSING.out.functional_annotation_folder
    motus_output = MOTUS_TAXONOMY.out.motus
    optional_tax_file_flag = NO_TAX_FILE_FLAG.out.created_file
    rna_count = RNA_PREDICTION.out.LSU_SSU_count
    sequence_categorisation_folder = MOVE_TO_SEQ_CAT_FOLDER.out.out
    stats = FUNCTIONAL_ANNOTATION_AND_POST_PROCESSING.out.stats
    taxonomy_summary_folder = RETURN_TAX_DIR.out.out

}
