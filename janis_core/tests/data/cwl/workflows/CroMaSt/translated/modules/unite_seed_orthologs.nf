nextflow.enable.dsl=2

process UNITE_SEED_ORTHOLOGS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/functional_annotation_and_post_processing/functional_annotation/eggnog/unite_seed_orthologs"
    memory "${params.after_qc.functional_annotation_and_post_processing.functional_annotation.eggnog.unite_seed_orthologs.memory}"

    input:
    path files
    val output_file_name

    output:
    stdout, emit: result

    script:
    def files = files.join(' ')
    """
    cat \
    ${files} \
    > stdout.txt \
    """

}
