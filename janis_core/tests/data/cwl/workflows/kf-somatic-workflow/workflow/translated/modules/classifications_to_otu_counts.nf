nextflow.enable.dsl=2

process CLASSIFICATIONS_TO_OTU_COUNTS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/rna_prediction/classify_ssus/classifications_to_otu_counts"

    input:
    path query, stageAs: 'query'
    val label
    val otu_table

    output:
    path "{inputs.query.basename}.tsv", emit: otu_tsv
    path "{inputs.query.basename}.notaxid.tsv", optional: true, emit: otu_tsv_notaxid
    path "{inputs.query.basename}.txt", emit: otu_txt

    script:
    def taxid_flag = params.after_qc.rna_prediction.classify_ssus.classifications_to_otu_counts_taxid_flag == false ? "" : "--taxid"
    """
    mapseq2biom.pl \
    --query ${query} \
    --label ${label} \
    --otuTable ${otu_table} \
    ${taxid_flag} \
    --krona "${query.name}.txt" \
    --notaxidfile "${query.name}.notaxid.tsv" \
    --outfile "${query.name}.tsv" \
    """

}


process CLASSIFICATIONS_TO_OTU_COUNTS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/rna_prediction/classify_ssus/classifications_to_otu_counts"

    input:
    path query, stageAs: 'query'
    val label
    val otu_table

    output:
    path "{inputs.query.basename}.tsv", emit: otu_tsv
    path "{inputs.query.basename}.notaxid.tsv", optional: true, emit: otu_tsv_notaxid
    path "{inputs.query.basename}.txt", emit: otu_txt

    script:
    def taxid_flag = params.after_qc.rna_prediction.classify_ssus.classifications_to_otu_counts_taxid_flag == false ? "" : "--taxid"
    """
    mapseq2biom.pl \
    --query ${query} \
    --label ${label} \
    --otuTable ${otu_table} \
    ${taxid_flag} \
    --krona "${query.name}.txt" \
    --notaxidfile "${query.name}.notaxid.tsv" \
    --outfile "${query.name}.tsv" \
    """

}


process CLASSIFICATIONS_TO_OTU_COUNTS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/rna_prediction/classify_ssus/classifications_to_otu_counts"

    input:
    path query, stageAs: 'query'
    val label
    val otu_table

    output:
    path "{inputs.query.basename}.tsv", emit: otu_tsv
    path "{inputs.query.basename}.notaxid.tsv", optional: true, emit: otu_tsv_notaxid
    path "{inputs.query.basename}.txt", emit: otu_txt

    script:
    def taxid_flag = params.after_qc.rna_prediction.classify_ssus.classifications_to_otu_counts_taxid_flag == false ? "" : "--taxid"
    """
    mapseq2biom.pl \
    --query ${query} \
    --label ${label} \
    --otuTable ${otu_table} \
    ${taxid_flag} \
    --krona "${query.name}.txt" \
    --notaxidfile "${query.name}.notaxid.tsv" \
    --outfile "${query.name}.tsv" \
    """

}
