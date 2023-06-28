nextflow.enable.dsl=2

process NO_TAX_FILE_FLAG {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/no_tax_file_flag"

    output:
    path "inputs.filename", emit: created_file

    script:
    """
    touch \
    ${params.after_qc.no_tax_file_flag_filename} \
    """

}


process NO_TAX_FILE_FLAG {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/no_tax_file_flag"

    output:
    path "inputs.filename", emit: created_file

    script:
    """
    touch \
    ${params.after_qc.no_tax_file_flag_filename} \
    """

}


process NO_TAX_FILE_FLAG {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/no_tax_file_flag"

    output:
    path "inputs.filename", emit: created_file

    script:
    """
    touch \
    ${params.after_qc.no_tax_file_flag_filename} \
    """

}
