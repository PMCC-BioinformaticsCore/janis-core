nextflow.enable.dsl=2

process SPLIT_FASTA {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/rna_prediction/find_ribosomal_ncrnas/cmsearch_raw_data/split_fasta"
    memory "${params.after_qc.rna_prediction.find_ribosomal_nc_rnas.cmsearch_raw_data.split_fasta.memory}"

    input:
    path seqs, stageAs: 'seqs'

    output:
    path "{inputs.seqs.basename}.*", emit: chunks

    script:
    """
    esl-ssplit.sh \
    ${params.after_qc.rna_prediction.find_ribosomal_nc_rnas.cmsearch_raw_data.split_fasta_number_of_output_files} \
    ${params.after_qc.rna_prediction.find_ribosomal_nc_rnas.cmsearch_raw_data.split_fasta_same_number_of_residues} \
    ${seqs} \
    ${params.after_qc.rna_prediction.find_ribosomal_nc_rnas.cmsearch_raw_data.split_fasta_chunk_size} \
    """

}


process SPLIT_FASTA {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/rna_prediction/find_ribosomal_ncrnas/cmsearch_raw_data/split_fasta"
    memory "${params.after_qc.rna_prediction.find_ribosomal_nc_rnas.cmsearch_raw_data.split_fasta.memory}"

    input:
    path seqs, stageAs: 'seqs'

    output:
    path "{inputs.seqs.basename}.*", emit: chunks

    script:
    """
    esl-ssplit.sh \
    ${params.after_qc.rna_prediction.find_ribosomal_nc_rnas.cmsearch_raw_data.split_fasta_number_of_output_files} \
    ${params.after_qc.rna_prediction.find_ribosomal_nc_rnas.cmsearch_raw_data.split_fasta_same_number_of_residues} \
    ${seqs} \
    ${params.after_qc.rna_prediction.find_ribosomal_nc_rnas.cmsearch_raw_data.split_fasta_chunk_size} \
    """

}
