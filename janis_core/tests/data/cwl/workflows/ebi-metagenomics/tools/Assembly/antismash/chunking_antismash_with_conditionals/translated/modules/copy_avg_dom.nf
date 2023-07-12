nextflow.enable.dsl=2

process COPY_AVG_DOM {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/unmapped_from_pfam/copy_avg_dom"

    input:
    path avg_unp_dom

    output:
    path "inputs.dir_cp", emit: dir_unp_dom

    script:
    """
    python script.py \
    """

}
