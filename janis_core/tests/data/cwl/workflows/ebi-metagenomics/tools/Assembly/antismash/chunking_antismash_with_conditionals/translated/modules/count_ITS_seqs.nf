nextflow.enable.dsl=2

process COUNT_ITS_SEQS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/its/count_its_seqs"

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
