nextflow.enable.dsl=2

process COUNT_SSU_FASTA {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/rna_prediction/count_ssu_fasta"

    input:
    path sequences, stageAs: 'sequences'

    output:
    val "grepcount", emit: count

    script:
    """
    grep -c ^> \
    ${sequences} \
    | \
    cat \
    > grepcount \
    """

}


process COUNT_SSU_FASTA {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/rna_prediction/count_ssu_fasta"

    input:
    path sequences, stageAs: 'sequences'

    output:
    val "grepcount", emit: count

    script:
    """
    grep -c ^> \
    ${sequences} \
    | \
    cat \
    > grepcount \
    """

}


process COUNT_SSU_FASTA {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/rna_prediction/count_ssu_fasta"

    input:
    path sequences, stageAs: 'sequences'

    output:
    val "grepcount", emit: count

    script:
    """
    grep -c ^> \
    ${sequences} \
    | \
    cat \
    > grepcount \
    """

}
