nextflow.enable.dsl=2

process CHUNKING_FASTA {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/antismash/chunking/chunking_fasta"

    input:
    path seqs, stageAs: 'seqs'
    val chunk_size

    output:
    path "{inputs.seqs.basename}.*", emit: chunks

    script:
    """
    esl-ssplit.sh \
    ${params.after_qc.antismash.chunking.chunking_fasta_number_of_output_files} \
    ${params.after_qc.antismash.chunking.chunking_fasta_same_number_of_residues} \
    ${seqs} \
    ${chunk_size} \
    """

}
