nextflow.enable.dsl=2

process MAPSEQ {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/rna_prediction/classify_ssus/mapseq"

    input:
    tuple path(primary), path(mscluster)
    path prefix, stageAs: 'prefix'
    path sequences, stageAs: 'sequences'
    val taxonomy

    output:
    stdout, emit: classifications

    script:
    """
    mapseq \
    -nthreads \
    8 \
    -tophits \
    80 \
    -topotus \
    40 \
    -outfmt \
    simple \
    ${sequences} \
    ${primary} \
    ${taxonomy} \
    > "${${prefix}.baseName}_${primary.name}.mseq" \
    """

}


process MAPSEQ {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/rna_prediction/classify_ssus/mapseq"

    input:
    tuple path(primary), path(mscluster)
    path prefix, stageAs: 'prefix'
    path sequences, stageAs: 'sequences'
    val taxonomy

    output:
    stdout, emit: classifications

    script:
    """
    mapseq \
    -nthreads \
    8 \
    -tophits \
    80 \
    -topotus \
    40 \
    -outfmt \
    simple \
    ${sequences} \
    ${primary} \
    ${taxonomy} \
    > "${${prefix}.baseName}_${primary.name}.mseq" \
    """

}


process MAPSEQ {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/rna_prediction/classify_ssus/mapseq"

    input:
    tuple path(primary), path(mscluster)
    path prefix, stageAs: 'prefix'
    path sequences, stageAs: 'sequences'
    val taxonomy

    output:
    stdout, emit: classifications

    script:
    """
    mapseq \
    -nthreads \
    8 \
    -tophits \
    80 \
    -topotus \
    40 \
    -outfmt \
    simple \
    ${sequences} \
    ${primary} \
    ${taxonomy} \
    > "${${prefix}.baseName}_${primary.name}.mseq" \
    """

}
