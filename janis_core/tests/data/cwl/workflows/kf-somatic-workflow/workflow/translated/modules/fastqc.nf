nextflow.enable.dsl=2

process FASTQC {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/fastqc"

    input:
    path fastq

    output:
    path "FASTQC/*.html", emit: html_files
    path "FASTQC/*.zip", emit: zip_files

    script:
    def fastq = fastq ? fastq.join(' ') : ""
    """
    fastqc \
    --threads 1 \
    --outdir \
    FASTQC \
    ${fastq} \
    """

}


process FASTQC {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/fastqc"

    input:
    path fastq

    output:
    path "FASTQC/*.html", emit: html_files
    path "FASTQC/*.zip", emit: zip_files

    script:
    def fastq = fastq ? fastq.join(' ') : ""
    """
    fastqc \
    --threads 1 \
    --outdir \
    FASTQC \
    ${fastq} \
    """

}
