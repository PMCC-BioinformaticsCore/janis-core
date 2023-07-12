nextflow.enable.dsl=2

process SMETANA {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/smetana"

    input:
    path gem

    output:
    path "{inputs.identifier}_SMETANA*", emit: detailed_output_tsv

    script:
    def gem = gem.join(' ')
    def solver = params.solver ? "--solver ${params.solver}" : ""
    """
    bash script.sh \
    --output "${params.identifier}_SMETANA" \
    --verbose \
    --flavor fbc2 \
    ${solver} \
    ${gem} \
    """

}
