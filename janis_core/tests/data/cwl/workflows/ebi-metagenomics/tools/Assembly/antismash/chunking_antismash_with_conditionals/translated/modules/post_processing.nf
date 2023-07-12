nextflow.enable.dsl=2

process POST_PROCESSING {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/cgc/combined_gene_caller/post_processing"
    cpus "${params.after_qc.cgc.combined_gene_caller.post_processing.cpus}"
    memory "${params.after_qc.cgc.combined_gene_caller.post_processing.memory}"

    input:
    path masking_file, stageAs: 'masking_file'
    path predicted_proteins_fgs_faa, stageAs: 'predicted_proteins_fgs_faa'
    path predicted_proteins_fgs_ffn, stageAs: 'predicted_proteins_fgs_ffn'
    path predicted_proteins_fgs_out, stageAs: 'predicted_proteins_fgs_out'
    val basename

    output:
    path "{inputs.basename}.faa", emit: predicted_proteins
    path "{inputs.basename}.ffn", emit: predicted_seq
    path "stderr.txt", emit: stderr
    stdout, emit: stdout

    script:
    """
    unite_protein_predictions.py \
    --fgs-faa ${predicted_proteins_fgs_faa} \
    --fgs-ffn ${predicted_proteins_fgs_ffn} \
    --fgs-out ${predicted_proteins_fgs_out} \
    --mask ${masking_file} \
    --name ${basename} \
    2> stderr.txt \
    > stdout.txt \
    """

}
