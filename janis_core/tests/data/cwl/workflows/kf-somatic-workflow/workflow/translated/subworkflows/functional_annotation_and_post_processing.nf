nextflow.enable.dsl=2

include { CHANGE_FORMATS_AND_NAMES } from './change_formats_and_names'
include { DIAMOND } from './diamond'
include { FOLDER_FUNCTIONAL_ANNOTATION } from './folder_functional_annotation'
include { FUNCTIONAL_ANNOTATION } from './functional_annotation'
include { GENOME_PROPERTIES } from '../modules/genome_properties'
include { GFF } from '../modules/gff'
include { MOVE_ANTISMASH_SUMMARY_TO_PATHWAYS_SYSTEMS_FOLDER } from '../modules/move_antismash_summary_to_pathways_systems_folder'
include { MOVE_TO_PATHWAYS_SYSTEMS_FOLDER } from '../modules/move_to_pathways_systems_folder'
include { PATHWAYS } from './pathways'
include { CHUNKING_TSV } from '../modules/chunking_tsv'
include { GO_SUMMARY } from '../modules/go_summary'
include { HEADER_ADDITION } from '../modules/header_addition'
include { MOVE_TO_FUNCTIONAL_ANNOTATION_FOLDER } from '../modules/move_to_functional_annotation_folder'
include { PFAM } from '../modules/pfam'
include { WRITE_SUMMARIES } from './write_summaries'

workflow FUNCTIONAL_ANNOTATION_AND_POST_PROCESSING {

    take:
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
    ch_antismash_geneclusters_txt
    ch_cgc_results_faa
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
    ch_pathways_classes
    ch_pathways_names
    ch_protein_chunk_size_ips
    ch_protein_chunk_size_eggnog
    ch_protein_chunk_size_hmm
    ch_rna_prediction_ncrna

    main:
    CHANGE_FORMATS_AND_NAMES(
        ch_filtered_fasta,
        GENOME_PROPERTIES.out.summary,
        PATHWAYS.out.kegg_pathways_summary
    )

    DIAMOND(
        ch_uniref90_db_txt,
        ch_diamond_database_file,
        ch_filtered_fasta,
        ch_diamond_max_target_seqs,
        params.after_qc.functional_annotation_and_post_processing.diamond_output_format,
        ch_cgc_results_faa,
        params.after_qc.functional_annotation_and_post_processing.diamond_strand,
        params.after_qc.functional_annotation_and_post_processing.diamond_threads
    )

    FOLDER_FUNCTIONAL_ANNOTATION(
        FUNCTIONAL_ANNOTATION.out.ips_result,
        ch_cgc_results_faa,
        ch_diamond_header,
        DIAMOND.out.post_processing_output,
        FUNCTIONAL_ANNOTATION.out.eggnog_annotations,
        FUNCTIONAL_ANNOTATION.out.eggnog_orthologs,
        ch_filtered_fasta,
        ch_go_config,
        FUNCTIONAL_ANNOTATION.out.hmm_result,
        ch_hmmsearch_header,
        ch_ips_header,
        ch_ko_file,
        GFF.out.output_gff_gz,
        GFF.out.output_gff_index,
        ch_rna_prediction_ncrna,
        ch_antismash_geneclusters_txt
    )

    FUNCTIONAL_ANNOTATION(
        ch_cgc_results_faa,
        ch_egg_nog_data_dir,
        ch_egg_nog_db,
        ch_egg_nog_diamond_db,
        ch_hmm_name_database,
        ch_hmm_gathering_bit_score,
        ch_hmm_omit_alignment,
        ch_inter_pro_scan_applications,
        ch_inter_pro_scan_databases,
        ch_inter_pro_scan_output_format,
        ch_protein_chunk_size_ips,
        ch_protein_chunk_size_eggnog,
        ch_protein_chunk_size_hmm,
        ch_func_ann_names_hmmer,
        ch_func_ann_names_ips
    )

    GENOME_PROPERTIES(
        FUNCTIONAL_ANNOTATION.out.ips_result,
        ch_gp_flatfiles_path,
        ch_filtered_fasta
    )

    GFF(
        FUNCTIONAL_ANNOTATION.out.eggnog_annotations,
        ch_cgc_results_faa,
        FUNCTIONAL_ANNOTATION.out.ips_result,
        ch_filtered_fasta
    )

    MOVE_ANTISMASH_SUMMARY_TO_PATHWAYS_SYSTEMS_FOLDER(
        FOLDER_FUNCTIONAL_ANNOTATION.out.summary_antismash
    )

    MOVE_TO_PATHWAYS_SYSTEMS_FOLDER(
        PATHWAYS.out.kegg_contigs_summary.toList()
    )

    PATHWAYS(
        ch_filtered_fasta,
        ch_graphs,
        FUNCTIONAL_ANNOTATION.out.hmm_result,
        ch_filtered_fasta,
        ch_pathways_classes,
        ch_pathways_names
    )

    emit:
    functional_annotation_folder = FOLDER_FUNCTIONAL_ANNOTATION.out.functional_annotation_folder
    pathways_systems_folder = MOVE_TO_PATHWAYS_SYSTEMS_FOLDER.out.out
    pathways_systems_folder_antismash_summary = MOVE_ANTISMASH_SUMMARY_TO_PATHWAYS_SYSTEMS_FOLDER.out.summary_in_folder
    stats = FOLDER_FUNCTIONAL_ANNOTATION.out.stats

}


workflow FUNCTIONAL_ANNOTATION_AND_POST_PROCESSING {

    take:
    ch_hmm_database
    ch_hmm_gathering_bit_score
    ch_hmm_omit_alignment
    ch_inter_pro_scan_applications
    ch_inter_pro_scan_databases
    ch_inter_pro_scan_output_format
    ch_cgc_results_faa
    ch_filtered_fasta
    ch_func_ann_names_hmmer
    ch_func_ann_names_ips
    ch_go_config
    ch_hmmsearch_header
    ch_ips_header
    ch_ko_file
    ch_protein_chunk_size_ips
    ch_protein_chunk_size_hmm
    ch_rna_prediction_ncrna

    main:
    CHUNKING_TSV(
        HEADER_ADDITION.out.output_table.toList()
    )

    FUNCTIONAL_ANNOTATION(
        ch_cgc_results_faa,
        ch_hmm_database,
        ch_hmm_gathering_bit_score,
        ch_hmm_omit_alignment,
        ch_inter_pro_scan_applications,
        ch_inter_pro_scan_databases,
        ch_inter_pro_scan_output_format,
        ch_protein_chunk_size_ips,
        ch_protein_chunk_size_hmm,
        ch_func_ann_names_hmmer,
        ch_func_ann_names_ips
    )

    FUNCTIONAL_ANNOTATION(
        ch_cgc_results_faa,
        ch_hmm_database,
        ch_hmm_gathering_bit_score,
        ch_hmm_omit_alignment,
        ch_inter_pro_scan_applications,
        ch_inter_pro_scan_databases,
        ch_inter_pro_scan_output_format,
        ch_protein_chunk_size_ips,
        ch_protein_chunk_size_hmm,
        ch_func_ann_names_hmmer,
        ch_func_ann_names_ips
    )

    GO_SUMMARY(
        FUNCTIONAL_ANNOTATION.out.ips_result,
        ch_go_config,
        ch_filtered_fasta
    )

    HEADER_ADDITION(
        FUNCTIONAL_ANNOTATION.out.hmm_result,
        ch_hmmsearch_header
    )

    MOVE_TO_FUNCTIONAL_ANNOTATION_FOLDER(
        WRITE_SUMMARIES.out.summary_ips.toList()
    )

    PFAM(
        FUNCTIONAL_ANNOTATION.out.ips_result,
        ch_filtered_fasta
    )

    WRITE_SUMMARIES(
        ch_cgc_results_faa,
        FUNCTIONAL_ANNOTATION.out.hmm_result,
        FUNCTIONAL_ANNOTATION.out.ips_result,
        ch_ko_file,
        PFAM.out.annotations,
        ch_rna_prediction_ncrna
    )

    emit:
    functional_annotation_folder = MOVE_TO_FUNCTIONAL_ANNOTATION_FOLDER.out.out
    stats = WRITE_SUMMARIES.out.stats

}
