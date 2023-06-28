nextflow.enable.dsl=2

process CAT {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/its/cat"

    input:
    path files

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
