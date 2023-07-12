nextflow.enable.dsl=2

process CLEAN_HEADERS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/before_qc/clean_headers"

    input:
    path sequences, stageAs: 'sequences'

    output:
    stdout, emit: sequences_with_cleaned_headers

    script:
    """
    tr " /|<_;#" ------- \
    > "${${sequences}.baseName}.unfiltered_fasta" \
    < ${sequences} \
    """

}
