nextflow.enable.dsl=2

process KOFAMSCAN {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/kofamscan"

    input:
    path input, stageAs: 'input'

    output:
    path "{inputs.identifier}.kofamscan.ttl", emit: output

    script:
    def threads = params.threads ? params.threads : 3
    """
    java -Xmx5g -jar /SAPP-2.0.jar -kofamscan \
    -input ${input} \
    -threads ${threads} \
    -kolist /ko_list \
    -output "${params.identifier}.kofamscan.ttl" \
    -profile /profiles/prokaryote.hal \
    """

}
