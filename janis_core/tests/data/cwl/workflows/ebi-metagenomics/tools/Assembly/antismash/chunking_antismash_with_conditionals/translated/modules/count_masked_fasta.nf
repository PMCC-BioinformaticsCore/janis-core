nextflow.enable.dsl=2

process COUNT_MASKED_FASTA {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/its/count_masked_fasta"

    input:
    path sequences, stageAs: 'sequences'

    output:
    val "grepcount", emit: count

    script:
    """
    grep -c ^> \
    ${sequences} \
    | \
    cat \
    > grepcount \
    """

}
