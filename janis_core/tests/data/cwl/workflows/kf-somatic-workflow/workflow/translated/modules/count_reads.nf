nextflow.enable.dsl=2

process COUNT_READS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/before_qc/count_reads"

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
