nextflow.enable.dsl=2

process MOTUS_CLASSIFICATION {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/motus_taxonomy/motus_classification"
    cpus "${params.after_qc.motus_taxonomy.motus_classification.cpus}"
    memory "${params.after_qc.motus_taxonomy.motus_classification.memory}"

    input:
    path reads, stageAs: 'reads'

    output:
    stdout, emit: motu_taxonomy

    script:
    """
    motus \
    -t 4 \
    profile \
    -c \
    -q \
    -s ${reads} \
    > "${${reads}.baseName}.motus" \
    """

}
