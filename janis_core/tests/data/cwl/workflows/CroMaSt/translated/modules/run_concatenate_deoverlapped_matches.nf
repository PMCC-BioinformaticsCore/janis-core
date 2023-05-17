nextflow.enable.dsl=2

process RUN_CONCATENATE_DEOVERLAPPED_MATCHES {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/rna_prediction/find_ribosomal_ncrnas/cmsearch_raw_data/run_concatenate_deoverlapped_matches"
    memory "${params.after_qc.rna_prediction.find_ribosomal_nc_rnas.cmsearch_raw_data.run_concatenate_deoverlapped_matches.memory}"

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


process RUN_CONCATENATE_DEOVERLAPPED_MATCHES {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/rna_prediction/find_ribosomal_ncrnas/cmsearch_raw_data/run_concatenate_deoverlapped_matches"
    memory "${params.after_qc.rna_prediction.find_ribosomal_nc_rnas.cmsearch_raw_data.run_concatenate_deoverlapped_matches.memory}"

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


process RUN_CONCATENATE_DEOVERLAPPED_MATCHES {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/rna_prediction/find_ribosomal_ncrnas/cmsearch_raw_data/run_concatenate_deoverlapped_matches"
    memory "${params.after_qc.rna_prediction.find_ribosomal_nc_rnas.cmsearch_raw_data.run_concatenate_deoverlapped_matches.memory}"

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
