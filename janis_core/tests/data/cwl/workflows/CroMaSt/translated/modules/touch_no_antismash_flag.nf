nextflow.enable.dsl=2

process TOUCH_NO_ANTISMASH_FLAG {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/antismash/no_antismash_subwf/touch_no_antismash_flag"

    output:
    path "inputs.filename", emit: created_file

    script:
    """
    touch \
    ${params.after_qc.antismash.no_antismash_subwf.touch_no_antismash_flag_filename} \
    """

}
