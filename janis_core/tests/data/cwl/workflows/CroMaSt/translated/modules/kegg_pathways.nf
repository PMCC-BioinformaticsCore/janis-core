nextflow.enable.dsl=2

process KEGG_PATHWAYS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/functional_annotation_and_post_processing/pathways/kegg_pathways"
    memory "${params.after_qc.functional_annotation_and_post_processing.pathways.kegg_pathways.memory}"

    input:
    path input_table, stageAs: 'input_table'
    val graphs
    val outputname
    val pathways_classes
    val pathways_names

    output:
    stdout, emit: stdout
    path "*summary.kegg_contigs*", emit: summary_contigs
    path "*summary.kegg_pathways*", emit: summary_pathways

    script:
    """
    give_pathways.py \
    -i ${input_table} \
    -c ${pathways_classes} \
    -g ${graphs} \
    -n ${pathways_names} \
    -o ${outputname} \
    > stdout.txt \
    """

}
