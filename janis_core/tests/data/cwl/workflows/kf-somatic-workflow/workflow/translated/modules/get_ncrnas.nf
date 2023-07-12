nextflow.enable.dsl=2

process GET_NCRNAS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/other_ncrnas/get_ncrnas"

    input:
    tuple path(primary), path(ssi)
    path names_contain_subseq_coords, stageAs: 'names_contain_subseq_coords'

    output:
    stdout, emit: sequences

    script:
    """
    esl-sfetch \
    -Cf ${primary} \
    ${names_contain_subseq_coords} \
    > "${primary.name}_${names_contain_subseq_coords.name}.fasta" \
    """

}


process GET_NCRNAS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/other_ncrnas/get_ncrnas"

    input:
    tuple path(primary), path(ssi)
    path names_contain_subseq_coords, stageAs: 'names_contain_subseq_coords'

    output:
    stdout, emit: sequences

    script:
    """
    esl-sfetch \
    -Cf ${primary} \
    ${names_contain_subseq_coords} \
    > "${primary.name}_${names_contain_subseq_coords.name}.fasta" \
    """

}
