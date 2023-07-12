nextflow.enable.dsl=2

process HEADER_ADDITION {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/functional_annotation_and_post_processing/header_addition"

    input:
    path input_table, stageAs: 'input_table'
    val header

    output:
    stdout, emit: output_table

    script:
    """
    add_header \
    -i ${input_table} \
    -h ${header} \
    > ${${input_table}.baseName} \
    """

}
