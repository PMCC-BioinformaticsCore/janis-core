nextflow.enable.dsl=2

process WRITE_SUMMARIES {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/functional_annotation_and_post_processing/write_summaries/write_summaries"

    input:
    path ips_entry_maps, stageAs: 'ips_entry_maps'
    path ko_entry_maps, stageAs: 'ko_entry_maps'
    path pfam_entry_maps, stageAs: 'pfam_entry_maps'
    path antismash_entry_maps, stageAs: 'antismash_entry_maps'
    val ips_outname
    val ko_outname
    val pfam_outname
    val antismash_outname

    output:
    path "*summary.antismash", optional: true, emit: summary_antismash
    path "*summary.ips", emit: summary_ips
    path "*summary.ko", emit: summary_ko
    path "*summary.pfam", emit: summary_pfam

    script:
    def antismash_entry_maps = antismash_entry_maps ? "-a ${antismash_entry_maps}" : ""
    def antismash_outname = antismash_outname ? "--antismash-name ${antismash_outname}" : ""
    """
    write_summaries.py \
    ${antismash_entry_maps} \
    -i ${ips_entry_maps} \
    -k ${ko_entry_maps} \
    -p ${pfam_entry_maps} \
    ${antismash_outname} \
    --ips-name ${ips_outname} \
    --ko-name ${ko_outname} \
    --pfam-name ${pfam_outname} \
    """

}
