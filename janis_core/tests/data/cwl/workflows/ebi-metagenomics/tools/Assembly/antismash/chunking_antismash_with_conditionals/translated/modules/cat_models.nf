nextflow.enable.dsl=2

process CAT_MODELS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/rna_prediction/find_ribosomal_ncrnas/cmsearch_raw_data/cat_models"
    memory "${params.after_qc.rna_prediction.find_ribosomal_nc_rnas.cmsearch_raw_data.cat_models.memory}"

    input:
    path files

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


process CAT_MODELS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/rna_prediction/find_ribosomal_ncrnas/cmsearch_raw_data/cat_models"
    memory "${params.after_qc.rna_prediction.find_ribosomal_nc_rnas.cmsearch_raw_data.cat_models.memory}"

    input:
    path files

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
