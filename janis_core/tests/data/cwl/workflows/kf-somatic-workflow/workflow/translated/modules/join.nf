nextflow.enable.dsl=2

process JOIN {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/functional_annotation_and_post_processing/diamond/post_processing_uniref90/join"

    input:
    path input_diamond, stageAs: 'input_diamond'
    val filename
    val input_db

    output:
    stdout, emit: output_join

    script:
    """
    diamond_post_run_join.sh \
    -i ${input_diamond} \
    -d ${input_db} \
    > "${filename}_summary.diamond.without_header" \
    """

}
