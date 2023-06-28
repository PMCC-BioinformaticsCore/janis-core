nextflow.enable.dsl=2

process CLEAN_CLASSIFICATION {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/motus_taxonomy/clean_classification"

    input:
    path taxonomy, stageAs: 'taxonomy'

    output:
    path "*.tsv", emit: clean_annotations

    script:
    """
    clean_motus_output.sh \
    ${taxonomy} \
    """

}
