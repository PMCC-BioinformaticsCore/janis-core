nextflow.enable.dsl=2

process COUNT_CDS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/cgc/count_cds"

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
