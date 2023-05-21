nextflow.enable.dsl=2

process QC_FLAG {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/before_qc/qc_flag"
    memory "${params.before_qc.qc_flag.memory}"

    input:
    val qc_count

    output:
    path "QC-*", emit: qc_flag

    script:
    """
    touch \
    """

}


process QC_FLAG {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/before_qc/qc_flag"
    memory "${params.before_qc.qc_flag.memory}"

    input:
    val qc_count

    output:
    path "QC-*", emit: qc_flag

    script:
    """
    touch \
    """

}


process QC_FLAG {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/before_qc/qc_flag"
    memory "${params.before_qc.qc_flag.memory}"

    input:
    val qc_count

    output:
    path "QC-*", emit: qc_flag

    script:
    """
    touch \
    """

}
