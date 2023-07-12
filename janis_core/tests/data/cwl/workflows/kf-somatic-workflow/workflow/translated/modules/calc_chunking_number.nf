nextflow.enable.dsl=2

process CALC_CHUNKING_NUMBER {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/antismash/chunking/calc_chunking_number"

    input:
    path sequences, stageAs: 'sequences'
    val number

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
