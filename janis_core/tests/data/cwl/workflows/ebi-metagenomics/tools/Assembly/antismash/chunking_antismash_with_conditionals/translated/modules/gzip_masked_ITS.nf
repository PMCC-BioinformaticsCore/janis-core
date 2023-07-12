nextflow.enable.dsl=2

process GZIP_MASKED_ITS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/its/gzip_masked_its"
    cpus "${params.after_qc.its.gzip_masked_its.cpus}"
    memory "${params.after_qc.its.gzip_masked_its.memory}"

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
