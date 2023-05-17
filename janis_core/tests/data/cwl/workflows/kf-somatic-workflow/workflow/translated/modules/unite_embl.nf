nextflow.enable.dsl=2

process UNITE_EMBL {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/antismash/chunking/unite_embl"

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
