nextflow.enable.dsl=2

process TOUCH_EMPTY_FASTA {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/before_qc/trimming/touch_empty_fasta"

    output:
    path "inputs.filename", emit: created_file

    script:
    """
    touch \
    ${params.before_qc.trimming.touch_empty_fasta_filename} \
    """

}
