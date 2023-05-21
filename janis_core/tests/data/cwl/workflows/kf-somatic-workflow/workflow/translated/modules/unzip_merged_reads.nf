nextflow.enable.dsl=2

process UNZIP_MERGED_READS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/before_qc/overlap_reads/unzip_merged_reads"

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


process UNZIP_MERGED_READS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/before_qc/overlap_reads/unzip_merged_reads"

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
