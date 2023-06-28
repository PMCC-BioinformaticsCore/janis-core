nextflow.enable.dsl=2

include { COMBINE } from '../modules/combine'
include { COMBINED_GENE_CALLER } from './combined_gene_caller'
include { COUNT_CDS } from '../modules/count_cds'
include { SPLIT_SEQS } from '../modules/split_seqs'

workflow CGC {

    take:
    ch_chunk_size
    ch_input_fasta
    ch_maskfile
    ch_postfixes

    main:
    COMBINE(
        COMBINED_GENE_CALLER.out.predicted_proteins.toList(),
        ch_input_fasta,
        ch_postfixes.flatten().first()
    )

    COMBINED_GENE_CALLER(
        SPLIT_SEQS.out.chunks.flatten().first(),
        ch_maskfile
    )

    COUNT_CDS(
        COMBINE.out.result
    )

    SPLIT_SEQS(
        ch_input_fasta,
        ch_chunk_size
    )

    emit:
    count_faa = COUNT_CDS.out.count
    results = COMBINE.out.result

}


workflow CGC {

    take:
    ch_chunk_size
    ch_input_fasta
    ch_maskfile
    ch_postfixes

    main:
    COMBINE(
        COMBINED_GENE_CALLER.out.predicted_proteins.toList(),
        ch_input_fasta,
        ch_postfixes.flatten().first()
    )

    COMBINE(
        COMBINED_GENE_CALLER.out.predicted_proteins.toList(),
        ch_input_fasta,
        ch_postfixes.flatten().first()
    )

    COMBINED_GENE_CALLER(
        SPLIT_SEQS.out.chunks.flatten().first(),
        ch_maskfile
    )

    COMBINED_GENE_CALLER(
        SPLIT_SEQS.out.chunks.flatten().first(),
        ch_maskfile
    )

    COUNT_CDS(
        COMBINE.out.result
    )

    COUNT_CDS(
        COMBINE.out.result
    )

    SPLIT_SEQS(
        ch_input_fasta,
        ch_chunk_size
    )

    SPLIT_SEQS(
        ch_input_fasta,
        ch_chunk_size
    )

    emit:
    count_faa = COUNT_CDS.out.count
    results = COMBINE.out.result

}
