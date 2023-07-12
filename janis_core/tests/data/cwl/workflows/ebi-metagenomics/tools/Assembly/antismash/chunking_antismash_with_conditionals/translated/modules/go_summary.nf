nextflow.enable.dsl=2

process GO_SUMMARY {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/functional_annotation_and_post_processing/go_summary"
    memory "${params.after_qc.functional_annotation_and_post_processing.go_summary.memory}"

    input:
    path inter_pro_scan_results, stageAs: 'inter_pro_scan_results'
    val config
    val output_name

    output:
    path "*.go", emit: go_summary
    path "*.go_slim", emit: go_summary_slim
    path "stderr.txt", emit: stderr
    stdout, emit: stdout

    script:
    """
    go_summary_pipeline-1.0.py \
    --input-file ${inter_pro_scan_results} \
    --config ${config} \
    --output-file ${output_name} \
    2> stderr.txt \
    > stdout.txt \
    """

}
