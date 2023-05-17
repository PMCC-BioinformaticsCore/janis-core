nextflow.enable.dsl=2

include { CHINKING_SC_FASTA_NUCLEOTIDE } from '../modules/chinking_SC_fasta_nucleotide'
include { CHINKING_FASTA_NUCLEOTIDE } from '../modules/chinking_fasta_nucleotide'
include { CHINKING_FASTA_PROTEINS } from '../modules/chinking_fasta_proteins'

workflow CHUNKING_FINAL {

    take:
    ch_faa
    ch_fasta
    ch_ffn
    ch_lsu
    ch_ssu

    main:
    CHINKING_SC_FASTA_NUCLEOTIDE(
        ch_lsu.toList()
    )

    CHINKING_FASTA_NUCLEOTIDE(
        ch_fasta.toList()
    )

    CHINKING_FASTA_PROTEINS(
        ch_faa.toList()
    )

    emit:
    SC_fasta_chunks = CHINKING_SC_FASTA_NUCLEOTIDE.out.chunks
    nucleotide_fasta_chunks = CHINKING_FASTA_NUCLEOTIDE.out.chunks
    protein_fasta_chunks = CHINKING_FASTA_PROTEINS.out.chunks

}


workflow CHUNKING_FINAL {

    take:
    ch_faa
    ch_fasta
    ch_ffn
    ch_lsu
    ch_ssu

    main:
    CHINKING_SC_FASTA_NUCLEOTIDE(
        ch_lsu.toList()
    )

    CHINKING_SC_FASTA_NUCLEOTIDE(
        ch_lsu.toList()
    )

    CHINKING_FASTA_NUCLEOTIDE(
        ch_fasta.toList()
    )

    CHINKING_FASTA_NUCLEOTIDE(
        ch_fasta.toList()
    )

    CHINKING_FASTA_PROTEINS(
        ch_faa.toList()
    )

    CHINKING_FASTA_PROTEINS(
        ch_faa.toList()
    )

    emit:
    SC_fasta_chunks = CHINKING_SC_FASTA_NUCLEOTIDE.out.chunks
    nucleotide_fasta_chunks = CHINKING_FASTA_NUCLEOTIDE.out.chunks
    protein_fasta_chunks = CHINKING_FASTA_PROTEINS.out.chunks

}
