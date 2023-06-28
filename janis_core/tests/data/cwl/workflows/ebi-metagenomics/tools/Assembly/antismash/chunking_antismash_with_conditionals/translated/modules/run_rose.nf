nextflow.enable.dsl=2

process RUN_ROSE {
    debug true
    container "biowardrobe2/rose:v0.0.2"
    publishDir "${params.outdir}/run_rose"

    input:
    tuple path(primary), path(bai)
    path annotation_file, stageAs: 'annotation_file'
    path binding_sites_file, stageAs: 'binding_sites_file'

    output:
    path "*Gateway_SuperEnhancers.bed", emit: gateway_super_enhancers_bed
    path "*Plot_points.png", emit: plot_points_pic

    script:
    """
    ROSE_main -o ./ \
    -i ${binding_sites_file} \
    -r ${primary} \
    -g ${annotation_file} \
    -s ${params.stitch_distance} \
    -t ${params.tss_distance} \
    """

}
