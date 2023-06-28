nextflow.enable.dsl=2

process FUNCTIONAL_STATS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/functional_annotation_and_post_processing/write_summaries/functional_stats"
    memory "${params.after_qc.functional_annotation_and_post_processing.write_summaries.functional_stats.memory}"

    input:
    path cds_file, stageAs: 'cds_file'
    path cmsearch_file, stageAs: 'cmsearch_file'
    path hmmscan, stageAs: 'hmmscan'
    path interproscan, stageAs: 'interproscan'
    path pfam, stageAs: 'pfam'
    val ko_file

    output:
    path "antismash*.yaml", optional: true, emit: antismash_yaml
    path "InterProScan*.yaml", emit: ips_yaml
    path "KO*.yaml", emit: ko_yaml
    path "pfam*.yaml", emit: pfam_yaml
    path "functional-annotation", emit: stats

    script:
    def antismash_file = params.after_qc.functional_annotation_and_post_processing.write_summaries.antismash_gene_clusters ? "-a ${params.after_qc.functional_annotation_and_post_processing.write_summaries.antismash_gene_clusters}" : ""
    """
    functional_stats.py \
    ${antismash_file} \
    -c ${cds_file} \
    -i ${interproscan} \
    -k ${hmmscan} \
    -p ${pfam} \
    -r ${cmsearch_file} \
    -ko ${ko_file} \
    """

}
