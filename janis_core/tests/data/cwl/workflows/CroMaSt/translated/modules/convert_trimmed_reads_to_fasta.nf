nextflow.enable.dsl=2

process CONVERT_TRIMMED_READS_TO_FASTA {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/before_qc/convert_trimmed_reads_to_fasta"
    memory "${params.before_qc.convert_trimmed_reads_to_fasta.memory}"

    input:
    path fastq, stageAs: 'fastq'

    output:
    path "*.unclean", emit: fasta

    script:
    """
    fastq_to_fasta.py \
    -i ${fastq} \
    -o "${${fastq}.baseName}.unclean" \
    """

}
