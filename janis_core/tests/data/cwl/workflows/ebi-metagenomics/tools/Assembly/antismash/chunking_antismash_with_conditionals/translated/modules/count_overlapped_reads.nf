nextflow.enable.dsl=2

process COUNT_OVERLAPPED_READS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/before_qc/trimming/count_overlapped_reads"
    cpus "${params.before_qc.trimming.count_overlapped_reads.cpus}"
    memory "${params.before_qc.trimming.count_overlapped_reads.memory}"

    input:
    path sequences, stageAs: 'sequences'

    output:
    val "data.txt", emit: count

    script:
    """
    count_lines.py \
    -f ${sequences} \
    -n ${params.before_qc.trimming.count_overlapped_reads_number} \
    """

}
