nextflow.enable.dsl=2

process GZIPPED_EMBL {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/antismash/chunking/gzipped_embl"
    cpus "${params.after_qc.antismash.chunking.gzipped_embl.cpus}"
    memory "${params.after_qc.antismash.chunking.gzipped_embl.memory}"

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
