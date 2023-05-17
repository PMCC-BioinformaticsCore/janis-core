nextflow.enable.dsl=2

process PFAM {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/functional_annotation_and_post_processing/pfam"
    memory "${params.after_qc.functional_annotation_and_post_processing.pfam.memory}"

    input:
    path interpro, stageAs: 'interpro'
    val outputname

    output:
    stdout, emit: annotations

    script:
    """
    \
    awk \
    /Pfam/ \
    ${interpro} \
    > ${outputname} \
    """

}
