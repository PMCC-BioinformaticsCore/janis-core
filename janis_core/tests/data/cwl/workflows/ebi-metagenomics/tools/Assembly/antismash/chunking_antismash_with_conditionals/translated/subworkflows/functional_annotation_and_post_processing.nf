nextflow.enable.dsl=2

include { CHUNKING_TSV } from '../modules/chunking_tsv'
include { FUNCTIONAL_ANNOTATION } from './functional_annotation'
include { GO_SUMMARY } from '../modules/go_summary'
include { HEADER_ADDITION } from '../modules/header_addition'
include { MOVE_TO_FUNCTIONAL_ANNOTATION_FOLDER } from '../modules/move_to_functional_annotation_folder'
include { PFAM } from '../modules/pfam'
include { WRITE_SUMMARIES } from './write_summaries'

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
