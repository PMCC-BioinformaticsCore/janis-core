nextflow.enable.dsl=2

include { CHUNKING_TSV } from '../modules/chunking_tsv'
include { COMPRESSION_FUNC_ANN } from '../modules/compression_func_ann'
include { GO_SUMMARY } from '../modules/go_summary'
include { HEADER_ADDITION } from '../modules/header_addition'
include { MOVE_TO_FUNCTIONAL_ANNOTATION_FOLDER } from '../modules/move_to_functional_annotation_folder'
include { PFAM } from '../modules/pfam'
include { WRITE_SUMMARIES } from './write_summaries'

workflow FOLDER_FUNCTIONAL_ANNOTATION {

    take:
    ch_ips_table
    ch_cds
    ch_diamond_header
    ch_diamond_table
    ch_eggnog_annotations
    ch_eggnog_orthologs
    ch_fasta
    ch_go_config
    ch_hmmscan_table
    ch_hmmsearch_header
    ch_ips_header
    ch_ko_file
    ch_output_gff_gz
    ch_output_gff_index
    ch_rna
    ch_antismash_geneclusters_txt

    main:
    CHUNKING_TSV(
        HEADER_ADDITION.out.output_table.toList()
    )

    COMPRESSION_FUNC_ANN(
        ch_eggnog_annotations
    )

    GO_SUMMARY(
        ch_ips_table,
        ch_go_config,
        ch_fasta
    )

    HEADER_ADDITION(
        ch_diamond_table,
        ch_diamond_header
    )

    MOVE_TO_FUNCTIONAL_ANNOTATION_FOLDER(
        ch_output_gff_gz.toList()
    )

    PFAM(
        ch_ips_table,
        ch_fasta
    )

    WRITE_SUMMARIES(
        ch_cds,
        ch_hmmscan_table,
        ch_ips_table,
        ch_ko_file,
        PFAM.out.annotations,
        ch_rna,
        ch_antismash_geneclusters_txt
    )

    emit:
    functional_annotation_folder = MOVE_TO_FUNCTIONAL_ANNOTATION_FOLDER.out.out
    stats = WRITE_SUMMARIES.out.stats
    summary_antismash = WRITE_SUMMARIES.out.summary_antismash

}
