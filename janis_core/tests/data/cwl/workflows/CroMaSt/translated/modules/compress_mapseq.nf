nextflow.enable.dsl=2

process COMPRESS_MAPSEQ {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/rna_prediction/classify_ssus/compress_mapseq"
    cpus "${params.after_qc.rna_prediction.classify_ssus.compress_mapseq.cpus}"
    memory "${params.after_qc.rna_prediction.classify_ssus.compress_mapseq.memory}"

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


process COMPRESS_MAPSEQ {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/rna_prediction/classify_ssus/compress_mapseq"
    cpus "${params.after_qc.rna_prediction.classify_ssus.compress_mapseq.cpus}"
    memory "${params.after_qc.rna_prediction.classify_ssus.compress_mapseq.memory}"

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


process COMPRESS_MAPSEQ {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/rna_prediction/classify_ssus/compress_mapseq"
    cpus "${params.after_qc.rna_prediction.classify_ssus.compress_mapseq.cpus}"
    memory "${params.after_qc.rna_prediction.classify_ssus.compress_mapseq.memory}"

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
