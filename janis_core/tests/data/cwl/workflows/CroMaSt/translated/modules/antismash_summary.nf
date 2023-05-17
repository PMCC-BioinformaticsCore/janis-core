nextflow.enable.dsl=2

process ANTISMASH_SUMMARY {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/antismash/chunking/antismash_summary"
    memory "${params.after_qc.antismash.chunking.antismash_summary.memory}"

    input:
    path geneclusters, stageAs: 'geneclusters'
    val glossary

    output:
    path "geneclusters-summary.txt", emit: reformatted_clusters

    script:
    """
    reformat_antismash.py \
    -g ${glossary} \
    -a ${geneclusters} \
    """

}
