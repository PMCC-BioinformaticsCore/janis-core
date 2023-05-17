nextflow.enable.dsl=2

process EGGNOG_ANNOTATION {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/functional_annotation_and_post_processing/functional_annotation/eggnog/eggnog_annotation"
    cpus "${params.after_qc.functional_annotation_and_post_processing.functional_annotation.eggnog.eggnog_annotation.cpus}"
    memory "${params.after_qc.functional_annotation_and_post_processing.functional_annotation.eggnog.eggnog_annotation.memory}"

    input:
    path annotate_hits_table, stageAs: 'annotate_hits_table'
    val cpu
    val data_dir
    val output

    output:
    path "{inputs.output}*annotations*", optional: true, emit: output_annotations
    path "{inputs.output}*orthologs*", optional: true, emit: output_orthologs

    script:
    def annotate_hits_table = annotate_hits_table ? "--annotate_hits_table ${annotate_hits_table}" : ""
    def cpu = cpu ? "--cpu ${cpu}" : ""
    def data_dir = data_dir ? "--data_dir ${data_dir}" : ""
    def no_file_comments = params.after_qc.functional_annotation_and_post_processing.functional_annotation.eggnog.eggnog_annotation_no_file_comments == false ? "" : "--no_file_comments"
    def output = output ? "-o ${output}" : ""
    """
    emapper.py \
    ${annotate_hits_table} \
    ${cpu} \
    ${data_dir} \
    ${output} \
    ${no_file_comments} \
    """

}
