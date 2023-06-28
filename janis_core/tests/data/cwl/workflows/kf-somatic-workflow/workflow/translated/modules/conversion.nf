nextflow.enable.dsl=2

process CONVERSION {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/conversion"

    input:
    path embl, stageAs: 'embl'

    output:
    path "{inputs.identifier}.ttl", emit: output

    script:
    def codon = params.codon ? "-codon ${params.codon}" : ""
    def embl = embl ? "-input ${embl}" : ""
    """
    java -Xmx5g -jar /SAPP-2.0.jar \
    ${embl} \
    ${codon} \
    -id ${params.identifier} \
    -output "${params.identifier}.ttl" \
    """

}
