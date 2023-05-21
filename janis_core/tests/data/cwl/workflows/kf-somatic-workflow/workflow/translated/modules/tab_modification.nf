nextflow.enable.dsl=2

process TAB_MODIFICATION {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/functional_annotation_and_post_processing/pathways/tab_modification"

    input:
    path input_table, stageAs: 'input_table'

    output:
    stdout, emit: output_with_tabs

    script:
    """
    sed /^#/d; s/ \+/\t/g \
    ${input_table} \
    > "${${input_table}.baseName}_tab.tbl" \
    """

}
