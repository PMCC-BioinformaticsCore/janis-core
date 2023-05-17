nextflow.enable.dsl=2

process GET_COORDS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/other_ncrnas/get_coords"
    memory "${params.after_qc.other_ncrnas.get_coords.memory}"

    input:
    path hits, stageAs: 'hits'
    val model

    output:
    path "*.RF*", emit: matches

    script:
    def model = model.join(' ')
    """
    pull_ncrnas.sh \
    ${hits} \
    ${model} \
    """

}


process GET_COORDS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/other_ncrnas/get_coords"
    memory "${params.after_qc.other_ncrnas.get_coords.memory}"

    input:
    path hits, stageAs: 'hits'
    val model

    output:
    path "*.RF*", emit: matches

    script:
    def model = model.join(' ')
    """
    pull_ncrnas.sh \
    ${hits} \
    ${model} \
    """

}
