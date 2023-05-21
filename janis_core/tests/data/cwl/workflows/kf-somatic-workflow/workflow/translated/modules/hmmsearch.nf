nextflow.enable.dsl=2

process HMMSEARCH {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/functional_annotation_and_post_processing/functional_annotation/run_hmmer/hmmsearch"

    input:
    path seqfile, stageAs: 'seqfile'
    val path_database
    val gathering_bit_score
    val omit_alignment

    output:
    path "*_hmmsearch.tbl", emit: output_table

    script:
    def gathering_bit_score = gathering_bit_score == false ? "" : "--cut_ga"
    def omit_alignment = omit_alignment == false ? "" : "--noali"
    """
    hmmsearch \
    --cpu 4 \
    -o /dev/null \
    ${omit_alignment} \
    --domtblout "${${seqfile}.baseName}_hmmsearch.tbl" \
    ${gathering_bit_score} \
    ${path_database} \
    ${seqfile} \
    > /dev/null \
    2> /dev/null \
    """

}


process HMMSEARCH {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/functional_annotation_and_post_processing/functional_annotation/run_hmmer/hmmsearch"

    input:
    path seqfile, stageAs: 'seqfile'
    val path_database
    val gathering_bit_score
    val omit_alignment

    output:
    path "*_hmmsearch.tbl", emit: output_table

    script:
    def gathering_bit_score = gathering_bit_score == false ? "" : "--cut_ga"
    def omit_alignment = omit_alignment == false ? "" : "--noali"
    """
    hmmsearch \
    --cpu 4 \
    -o /dev/null \
    ${omit_alignment} \
    --domtblout "${${seqfile}.baseName}_hmmsearch.tbl" \
    ${gathering_bit_score} \
    ${path_database} \
    ${seqfile} \
    > /dev/null \
    2> /dev/null \
    """

}
