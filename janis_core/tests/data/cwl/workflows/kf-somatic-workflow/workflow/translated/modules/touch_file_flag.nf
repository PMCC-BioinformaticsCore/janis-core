nextflow.enable.dsl=2

process TOUCH_FILE_FLAG {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/touch_file_flag"

    output:
    path "inputs.filename", emit: created_file

    script:
    """
    touch \
    ${params.touch_file_flag_filename} \
    """

}


process TOUCH_FILE_FLAG {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/touch_file_flag"

    output:
    path "inputs.filename", emit: created_file

    script:
    """
    touch \
    ${params.touch_file_flag_filename} \
    """

}


process TOUCH_FILE_FLAG {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/touch_file_flag"

    output:
    path "inputs.filename", emit: created_file

    script:
    """
    touch \
    ${params.touch_file_flag_filename} \
    """

}
