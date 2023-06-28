nextflow.enable.dsl=2

process VISUALIZE_OTU_COUNTS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/rna_prediction/classify_ssus/visualize_otu_counts"
    cpus "${params.after_qc.rna_prediction.classify_ssus.visualize_otu_counts.cpus}"
    memory "${params.after_qc.rna_prediction.classify_ssus.visualize_otu_counts.memory}"

    input:
    path otu_counts, stageAs: 'otu_counts'

    output:
    path "*.html", emit: otu_visualization

    script:
    """
    ktImportText \
    -o krona.html \
    ${otu_counts} \
    """

}


process VISUALIZE_OTU_COUNTS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/rna_prediction/classify_ssus/visualize_otu_counts"
    cpus "${params.after_qc.rna_prediction.classify_ssus.visualize_otu_counts.cpus}"
    memory "${params.after_qc.rna_prediction.classify_ssus.visualize_otu_counts.memory}"

    input:
    path otu_counts, stageAs: 'otu_counts'

    output:
    path "*.html", emit: otu_visualization

    script:
    """
    ktImportText \
    -o krona.html \
    ${otu_counts} \
    """

}
