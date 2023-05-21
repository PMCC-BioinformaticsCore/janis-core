nextflow.enable.dsl=2

process COMPRESSION_FUNC_ANN {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/functional_annotation_and_post_processing/folder_functional_annotation/compression_func_ann"

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
