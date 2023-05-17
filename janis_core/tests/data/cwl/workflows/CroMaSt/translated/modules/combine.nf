nextflow.enable.dsl=2

process COMBINE {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/functional_annotation_and_post_processing/functional_annotation/run_hmmer/combine"
    memory "${params.after_qc.functional_annotation_and_post_processing.functional_annotation.run_hmmer.combine.memory}"

    input:
    path files
    val output_file_name
    val postfix

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


process COMBINE {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/functional_annotation_and_post_processing/functional_annotation/run_hmmer/combine"
    memory "${params.after_qc.functional_annotation_and_post_processing.functional_annotation.run_hmmer.combine.memory}"

    input:
    path files
    val output_file_name
    val postfix

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
