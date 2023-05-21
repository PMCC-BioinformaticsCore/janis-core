nextflow.enable.dsl=2

process INTERPROSCAN {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/interproscan"
    memory "${params.interproscan.memory}"

    input:
    path input, stageAs: 'input'

    output:
    path "{inputs.identifier}.hdt", emit: output

    script:
    def cpu = params.threads ? params.threads : 2
    """
    java -Xmx5g -jar /SAPP-2.0.jar -interpro \
    -input ${input} \
    -cpu ${cpu} \
    -path ${params.interpro} \
    -a Pfam \
    -output "${params.identifier}.hdt" \
    """

}
