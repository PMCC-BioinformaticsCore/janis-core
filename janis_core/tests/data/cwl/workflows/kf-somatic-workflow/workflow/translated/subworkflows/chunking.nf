nextflow.enable.dsl=2

include { ANTISMASH_GFF } from '../modules/antismash_gff'
include { ANTISMASH_SUMMARY } from '../modules/antismash_summary'
include { CALC_CHUNKING_NUMBER } from '../modules/calc_chunking_number'
include { CHECK_VALUE } from '../modules/check_value'
include { CHUNKING_FASTA } from '../modules/chunking_fasta'
include { GZIPPED_EMBL } from '../modules/gzipped_embl'
include { GZIPPED_GBK } from '../modules/gzipped_gbk'
include { RENAME_CONTIGS } from '../modules/rename_contigs'
include { RENAME_GENECLUSTERS } from '../modules/rename_geneclusters'
include { RETURN_ANTISMASH_IN_FOLDER } from '../modules/return_antismash_in_folder'
include { RUN_ANTISMASH } from './run_antismash'
include { UNITE_EMBL } from '../modules/unite_embl'
include { UNITE_GBK } from '../modules/unite_gbk'
include { UNITE_GENECLUSTERS_TXT } from '../modules/unite_geneclusters_txt'

workflow CHUNKING {

    take:
    ch_clusters_glossary
    ch_filtered_fasta
    ch_final_folder_name
    ch_split_size

    main:
    ANTISMASH_GFF(
        UNITE_EMBL.out.result,
        ANTISMASH_SUMMARY.out.reformatted_clusters,
        ch_filtered_fasta
    )

    ANTISMASH_SUMMARY(
        UNITE_GENECLUSTERS_TXT.out.result,
        ch_clusters_glossary
    )

    CALC_CHUNKING_NUMBER(
        ch_filtered_fasta,
        ch_split_size
    )

    CHECK_VALUE(
        CALC_CHUNKING_NUMBER.out.count
    )

    CHUNKING_FASTA(
        ch_filtered_fasta,
        CHECK_VALUE.out.out
    )

    GZIPPED_EMBL(
        UNITE_EMBL.out.result
    )

    GZIPPED_GBK(
        UNITE_GBK.out.result
    )

    RENAME_CONTIGS(
        CHUNKING_FASTA.out.chunks.flatten().first(),
        ch_filtered_fasta,
        ch_filtered_fasta
    )

    RENAME_GENECLUSTERS(
        ANTISMASH_SUMMARY.out.reformatted_clusters,
        ch_filtered_fasta
    )

    RETURN_ANTISMASH_IN_FOLDER(
        ANTISMASH_GFF.out.output_gff_bgz.toList(),
        ch_final_folder_name
    )

    RUN_ANTISMASH(
        ch_filtered_fasta,
        RENAME_CONTIGS.out.renamed_contigs_in_chunks,
        RENAME_CONTIGS.out.names_table
    )

    UNITE_EMBL(
        RUN_ANTISMASH.out.antismash_embl.toList(),
        ch_filtered_fasta
    )

    UNITE_GBK(
        RUN_ANTISMASH.out.antismash_gbk.toList(),
        ch_filtered_fasta
    )

    UNITE_GENECLUSTERS_TXT(
        RUN_ANTISMASH.out.antismash_txt.toList()
    )

    emit:
    antismash_clusters = RENAME_GENECLUSTERS.out.renamed_file
    antismash_folder_chunking = RETURN_ANTISMASH_IN_FOLDER.out.out

}
