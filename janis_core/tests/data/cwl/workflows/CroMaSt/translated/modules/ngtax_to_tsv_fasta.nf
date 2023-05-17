nextflow.enable.dsl=2

process NGTAX_TO_TSV_FASTA {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/ngtax_to_tsv_fasta"

    input:
    path input, stageAs: 'input'
    path metadata, stageAs: 'metadata'

    output:
    path "*_asv.tsv", emit: physeq_asv
    path "*_met.tsv", emit: physeq_met
    path "*_seq.tsv", emit: physeq_seq
    path "*_tax.tsv", emit: physeq_tax
    path "*.picrust.fasta", emit: picrust_fasta
    path "*.picrust.tsv", emit: picrust_tsv

    script:
    """
    python3 /scripts/ngtax_to_tsv-fasta.py \
    -t \
    ${input} \
    -i \
    ${params.sample} \
    -f \
    ${params.fragment} \
    -m \
    ${metadata} \
    """

}
