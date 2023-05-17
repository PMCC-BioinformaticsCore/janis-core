nextflow.enable.dsl=2

process TOHDT {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/tohdt"

    input:
    path input, stageAs: 'input'

    output:
    path "inputs.output", emit: output

    script:
    """
    java -Xmx5g -jar /SAPP-2.0.jar -convert \
    -i ${input} \
    -o ${params.output} \
    """

}
