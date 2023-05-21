nextflow.enable.dsl=2

process SORTING {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/functional_annotation_and_post_processing/diamond/post_processing_uniref90/sorting"

    input:
    path input_table, stageAs: 'input_table'

    output:
    stdout, emit: output_sorted

    script:
    """
    sort -k2,2 \
    ${input_table} \
    > "${${input_table}.baseName}.sorted" \
    """

}
