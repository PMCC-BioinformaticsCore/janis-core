nextflow.enable.dsl=2

process PRODIGAL {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/prodigal"

    input:
    path input, stageAs: 'input'

    output:
    path "{inputs.identifier}.prodigal.ttl", emit: output

    script:
    """
    java -Xmx5g -jar /SAPP-2.0.jar -prodigal \
    -input ${input} \
    -output "${params.identifier}.prodigal.ttl" \
    """

}
