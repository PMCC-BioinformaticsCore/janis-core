nextflow.enable.dsl=2

process COMBINE_TRIMMED {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/before_qc/trimming/trimming/combine_trimmed"
    memory "${params.before_qc.trimming.trimming.combine_trimmed.memory}"

    input:
    path files
    val output_file_name

    output:
    stdout, emit: result

    script:
    def files = files.join(' ')
    """
    cat \
    ${files} \
    > stdout.txt \
    """

}
