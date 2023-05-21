nextflow.enable.dsl=2

process CREATE_CSV_GP {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/functional_annotation_and_post_processing/change_formats_and_names/create_csv_gp"
    memory "${params.after_qc.functional_annotation_and_post_processing.change_formats_and_names.create_csv_gp.memory}"

    input:
    path tab_sep_table, stageAs: 'tab_sep_table'
    val output_name

    output:
    path "inputs.output_name", emit: csv_result

    script:
    """
    make_csv.py \
    -i ${tab_sep_table} \
    -o ${output_name} \
    """

}
