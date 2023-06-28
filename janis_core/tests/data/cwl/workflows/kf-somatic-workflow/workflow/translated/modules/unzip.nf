nextflow.enable.dsl=2

process UNZIP {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/before_qc/unzip"

    input:
    path target_reads, stageAs: 'target_reads'

    output:
    stdout, emit: unzipped_file

    script:
    """
    gunzip -c \
    ${target_reads} \
    > stdout.txt \
    """

}
