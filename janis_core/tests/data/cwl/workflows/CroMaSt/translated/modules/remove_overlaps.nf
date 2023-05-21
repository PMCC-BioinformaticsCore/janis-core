nextflow.enable.dsl=2

process REMOVE_OVERLAPS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/rna_prediction/find_ribosomal_ncrnas/cmsearch_raw_data/remove_overlaps"
    cpus "${params.after_qc.rna_prediction.find_ribosomal_nc_rnas.cmsearch_raw_data.remove_overlaps.cpus}"
    memory "${params.after_qc.rna_prediction.find_ribosomal_nc_rnas.cmsearch_raw_data.remove_overlaps.memory}"

    input:
    path cmsearch_matches, stageAs: 'cmsearch_matches'
    val clan_information

    output:
    path "*.deoverlapped", emit: deoverlapped_matches

    script:
    def clan_information = clan_information ? "--clanin ${clan_information}" : ""
    """
    cmsearch-deoverlap.pl \
    ${clan_information} \
    ${cmsearch_matches} \
    """

}


process REMOVE_OVERLAPS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/rna_prediction/find_ribosomal_ncrnas/cmsearch_raw_data/remove_overlaps"
    cpus "${params.after_qc.rna_prediction.find_ribosomal_nc_rnas.cmsearch_raw_data.remove_overlaps.cpus}"
    memory "${params.after_qc.rna_prediction.find_ribosomal_nc_rnas.cmsearch_raw_data.remove_overlaps.memory}"

    input:
    path cmsearch_matches, stageAs: 'cmsearch_matches'
    val clan_information

    output:
    path "*.deoverlapped", emit: deoverlapped_matches

    script:
    def clan_information = clan_information ? "--clanin ${clan_information}" : ""
    """
    cmsearch-deoverlap.pl \
    ${clan_information} \
    ${cmsearch_matches} \
    """

}


process REMOVE_OVERLAPS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/rna_prediction/find_ribosomal_ncrnas/cmsearch_raw_data/remove_overlaps"
    cpus "${params.after_qc.rna_prediction.find_ribosomal_nc_rnas.cmsearch_raw_data.remove_overlaps.cpus}"
    memory "${params.after_qc.rna_prediction.find_ribosomal_nc_rnas.cmsearch_raw_data.remove_overlaps.memory}"

    input:
    path cmsearch_matches, stageAs: 'cmsearch_matches'
    val clan_information

    output:
    path "*.deoverlapped", emit: deoverlapped_matches

    script:
    def clan_information = clan_information ? "--clanin ${clan_information}" : ""
    """
    cmsearch-deoverlap.pl \
    ${clan_information} \
    ${cmsearch_matches} \
    """

}
