nextflow.enable.dsl=2

process MOVE_ANTISMASH_SUMMARY_TO_PATHWAYS_SYSTEMS_FOLDER {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/functional_annotation_and_post_processing/move_antismash_summary_to_pathways_systems_folder"

    input:
    path antismash_summary, stageAs: 'antismash_summary'

    output:
    path "inputs.folder_name", emit: summary_in_folder

    script:
    def antismash_summary = antismash_summary ? "-a ${antismash_summary}" : ""
    """
    move_antismash_summary.py \
    ${antismash_summary} \
    -f ${params.after_qc.functional_annotation_and_post_processing.move_antismash_summary_to_pathways_systems_folder_folder_name} \
    """

}
