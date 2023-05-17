nextflow.enable.dsl=2

process GENOME_PROPERTIES {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/functional_annotation_and_post_processing/genome_properties"
    memory "${params.after_qc.functional_annotation_and_post_processing.genome_properties.memory}"

    input:
    path input_tsv_file, stageAs: 'input_tsv_file'
    val flatfiles_path
    val name

    output:
    path "JSON*{inputs.name}", optional: true, emit: json
    path "stderr.txt", emit: stderr
    stdout, emit: stdout
    path "SUMMARY*{inputs.name}", emit: summary
    path "TABLE*{inputs.name}", optional: true, emit: table

    script:
    def name = name ? "-name ${name}" : ""
    """
    assign_genome_properties.pl \
    -matches ${input_tsv_file} \
    -gpdir ${flatfiles_path} \
    -gpff ${params.after_qc.functional_annotation_and_post_processing.genome_properties_gp_txt} \
    ${name} \
    -all \
    -outfiles table \
    -outfiles web_json \
    -outfiles summary \
    2> stderr.txt \
    > stdout.txt \
    """

}
