nextflow.enable.dsl=2

process PARSING_HMMSCAN {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/functional_annotation_and_post_processing/pathways/parsing_hmmscan"

    input:
    path fasta, stageAs: 'fasta'
    path table, stageAs: 'table'

    output:
    path "*_parsed*", emit: output_table

    script:
    """
    parsing_hmmscan.py \
    -f ${fasta} \
    -i ${table} \
    2> stderr.txt \
    > stdout.txt \
    """

}
