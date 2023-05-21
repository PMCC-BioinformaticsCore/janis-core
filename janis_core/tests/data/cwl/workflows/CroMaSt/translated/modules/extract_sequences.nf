nextflow.enable.dsl=2

process EXTRACT_SEQUENCES {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/rna_prediction/extract_sequences"
    cpus "${params.after_qc.rna_prediction.extract_sequences.cpus}"
    memory "${params.after_qc.rna_prediction.extract_sequences.memory}"

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


process EXTRACT_SEQUENCES {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/rna_prediction/extract_sequences"
    cpus "${params.after_qc.rna_prediction.extract_sequences.cpus}"
    memory "${params.after_qc.rna_prediction.extract_sequences.memory}"

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


process EXTRACT_SEQUENCES {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/rna_prediction/extract_sequences"
    cpus "${params.after_qc.rna_prediction.extract_sequences.cpus}"
    memory "${params.after_qc.rna_prediction.extract_sequences.memory}"

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
