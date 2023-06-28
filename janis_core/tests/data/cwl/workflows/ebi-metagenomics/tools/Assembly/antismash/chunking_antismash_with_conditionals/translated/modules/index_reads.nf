nextflow.enable.dsl=2

process INDEX_READS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/other_ncrnas/index_reads"
    cpus "${params.after_qc.other_ncrnas.index_reads.cpus}"
    memory "${params.after_qc.other_ncrnas.index_reads.memory}"

    input:
    path sequences, stageAs: 'sequences'

    output:
    tuple path("folder/{inputs.sequences.basename}"), path("*.ssi"), emit: sequences_with_index

    script:
    """
    esl-index.sh \
    -f ${sequences} \
    """

}
