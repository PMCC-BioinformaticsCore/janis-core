nextflow.enable.dsl=2

process GZIPPED_GBK {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/antismash/chunking/gzipped_gbk"
    cpus "${params.after_qc.antismash.chunking.gzipped_gbk.cpus}"
    memory "${params.after_qc.antismash.chunking.gzipped_gbk.memory}"

    input:
    path uncompressed_file, stageAs: 'uncompressed_file'

    output:
    stdout, emit: compressed_file

    script:
    """
    pigz \
    -p \
    16 \
    -c \
    ${uncompressed_file} \
    > "${uncompressed_file.name}.gz" \
    """

}
