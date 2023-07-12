nextflow.enable.dsl=2

process RUN_CONCATENATE_MATCHES {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/rna_prediction/find_ribosomal_ncrnas/cmsearch_raw_data/run_concatenate_matches"

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


process RUN_CONCATENATE_MATCHES {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/rna_prediction/find_ribosomal_ncrnas/cmsearch_raw_data/run_concatenate_matches"

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


process RUN_CONCATENATE_MATCHES {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/rna_prediction/find_ribosomal_ncrnas/cmsearch_raw_data/run_concatenate_matches"

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
