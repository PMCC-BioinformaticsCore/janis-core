nextflow.enable.dsl=2

process COUNT_READS_AFTER_FILTERING {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/antismash/filtering/count_reads_after_filtering"

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
