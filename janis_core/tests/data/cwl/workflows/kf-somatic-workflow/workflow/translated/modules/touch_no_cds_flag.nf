nextflow.enable.dsl=2

process TOUCH_NO_CDS_FLAG {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/touch_no_cds_flag"

    output:
    path "inputs.filename", emit: created_file

    script:
    """
    touch \
    ${params.touch_no_cds_flag_filename} \
    """

}


process TOUCH_NO_CDS_FLAG {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/touch_no_cds_flag"

    output:
    path "inputs.filename", emit: created_file

    script:
    """
    touch \
    ${params.touch_no_cds_flag_filename} \
    """

}
