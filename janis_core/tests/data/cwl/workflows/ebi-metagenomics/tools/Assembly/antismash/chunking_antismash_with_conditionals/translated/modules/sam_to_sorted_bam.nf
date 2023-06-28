nextflow.enable.dsl=2

process SAM_TO_SORTED_BAM {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/sam_to_sorted_bam"

    input:
    path sam, stageAs: 'sam'

    output:
    path "{inputs.identifier}.sorted.bam", emit: sortedbam

    script:
    def threads = params.threads ? params.threads : 1
    """
    bash -x script.sh \
    ${params.identifier} \
    ${threads} \
    ${sam} \
    """

}
