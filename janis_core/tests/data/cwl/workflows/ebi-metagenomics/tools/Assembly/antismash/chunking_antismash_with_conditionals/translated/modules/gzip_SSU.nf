nextflow.enable.dsl=2

process GZIP_SSU {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/gzip_ssu"
    cpus "${params.after_qc.gzip_ssu.cpus}"
    memory "${params.after_qc.gzip_ssu.memory}"

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
