nextflow.enable.dsl=2

process COMPRESSION {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/compression"

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


process COMPRESSION {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/compression"

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
