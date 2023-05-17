nextflow.enable.dsl=2

process UNITE_GENECLUSTERS_TXT {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/antismash/chunking/unite_geneclusters_txt"
    memory "${params.after_qc.antismash.chunking.unite_geneclusters_txt.memory}"

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
