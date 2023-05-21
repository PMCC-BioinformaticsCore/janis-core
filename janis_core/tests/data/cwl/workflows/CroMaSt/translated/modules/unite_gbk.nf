nextflow.enable.dsl=2

process UNITE_GBK {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/antismash/chunking/unite_gbk"
    memory "${params.after_qc.antismash.chunking.unite_gbk.memory}"

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
