nextflow.enable.dsl=2

process SPLIT_SEQS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/functional_annotation_and_post_processing/functional_annotation/split_seqs"
    memory "${params.after_qc.functional_annotation_and_post_processing.functional_annotation.split_seqs.memory}"

    input:
    path seqs, stageAs: 'seqs'
    val chunk_size

    output:
    path "*_*", emit: chunks

    script:
    """
    split_to_chunks.py \
    -i ${seqs} \
    -s ${chunk_size} \
    > /dev/null \
    2> /dev/null \
    """

}
