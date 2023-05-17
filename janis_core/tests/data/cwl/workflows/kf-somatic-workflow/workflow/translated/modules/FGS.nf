nextflow.enable.dsl=2

process FGS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/cgc/combined_gene_caller/fgs"

    input:
    path input_fasta, stageAs: 'input_fasta'
    val output

    output:
    path "{inputs.output}.faa", emit: predicted_proteins_faa
    path "{inputs.output}.ffn", emit: predicted_proteins_ffn
    path "{inputs.output}.out", emit: predicted_proteins_out
    path "stderr.txt", emit: stderr
    stdout, emit: stdout

    script:
    """
    run_FGS.sh \
    -i ${input_fasta} \
    -o ${output} \
    -s ${params.after_qc.cgc.combined_gene_caller.fgs_seq_type} \
    2> stderr.txt \
    > stdout.txt \
    """

}


process FGS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/cgc/combined_gene_caller/fgs"

    input:
    path input_fasta, stageAs: 'input_fasta'
    val output

    output:
    path "{inputs.output}.faa", emit: predicted_proteins_faa
    path "{inputs.output}.ffn", emit: predicted_proteins_ffn
    path "{inputs.output}.out", emit: predicted_proteins_out
    path "stderr.txt", emit: stderr
    stdout, emit: stdout

    script:
    """
    run_FGS.sh \
    -i ${input_fasta} \
    -o ${output} \
    -s ${params.after_qc.cgc.combined_gene_caller.fgs_seq_type} \
    2> stderr.txt \
    > stdout.txt \
    """

}
