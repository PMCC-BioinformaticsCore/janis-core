nextflow.enable.dsl=2

process EGGNOG_HOMOLOGY_SEARCHES {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/functional_annotation_and_post_processing/functional_annotation/eggnog/eggnog_homology_searches"

    input:
    path fasta_file, stageAs: 'fasta_file'
    val cpu
    val data_dir
    val db
    val db_diamond
    val output

    output:
    path "{inputs.output}*annotations*", optional: true, emit: output_annotations
    path "{inputs.output}*orthologs*", optional: true, emit: output_orthologs

    script:
    def cpu = cpu ? "--cpu ${cpu}" : ""
    def data_dir = data_dir ? "--data_dir ${data_dir}" : ""
    def db = db ? "--database ${db}" : ""
    def db_diamond = db_diamond ? "--dmnd_db ${db_diamond}" : ""
    def fasta_file = fasta_file ? "-i ${fasta_file}" : ""
    def mode = params.after_qc.functional_annotation_and_post_processing.functional_annotation.eggnog.eggnog_homology_searches_mode ? "-m ${params.after_qc.functional_annotation_and_post_processing.functional_annotation.eggnog.eggnog_homology_searches_mode}" : ""
    def no_annot = params.after_qc.functional_annotation_and_post_processing.functional_annotation.eggnog.eggnog_homology_searches_no_annot == false ? "" : "--no_annot"
    def no_file_comments = params.after_qc.functional_annotation_and_post_processing.functional_annotation.eggnog.eggnog_homology_searches_no_file_comments == false ? "" : "--no_file_comments"
    def output = output ? "-o ${output}" : ""
    """
    emapper.py \
    ${fasta_file} \
    ${cpu} \
    ${data_dir} \
    ${db} \
    ${db_diamond} \
    ${mode} \
    ${output} \
    ${no_annot} \
    ${no_file_comments} \
    """

}
