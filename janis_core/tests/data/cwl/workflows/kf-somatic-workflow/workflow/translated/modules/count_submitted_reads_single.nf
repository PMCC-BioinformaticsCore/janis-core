nextflow.enable.dsl=2

process COUNT_SUBMITTED_READS_SINGLE {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/before_qc/overlap_reads/count_submitted_reads_single"

    input:
    path sequences, stageAs: 'sequences'

    output:
    val "data.txt", emit: count

    script:
    """
    count_lines.py \
    -f ${sequences} \
    -n ${params.before_qc.overlap_reads.count_submitted_reads_single_number} \
    """

}


process COUNT_SUBMITTED_READS_SINGLE {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/before_qc/overlap_reads/count_submitted_reads_single"

    input:
    path sequences, stageAs: 'sequences'

    output:
    val "data.txt", emit: count

    script:
    """
    count_lines.py \
    -f ${sequences} \
    -n ${params.before_qc.overlap_reads.count_submitted_reads_single_number} \
    """

}
