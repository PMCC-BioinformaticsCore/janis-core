nextflow.enable.dsl=2

process GZIP_FILES {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/other_ncrnas/gzip_files"
    cpus "${params.after_qc.other_ncrnas.gzip_files.cpus}"
    memory "${params.after_qc.other_ncrnas.gzip_files.memory}"

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
