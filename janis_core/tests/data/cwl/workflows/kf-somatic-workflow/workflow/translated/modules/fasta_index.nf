nextflow.enable.dsl=2

process FASTA_INDEX {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/fasta_index"

    input:
    path fasta, stageAs: 'fasta'

    output:
    path "{inputs.fasta.basename}.bgz.gzi", emit: bgz_index
    path "{inputs.fasta.basename}.bgz", emit: fasta_bgz
    path "{inputs.fasta.basename}.bgz.fai", emit: fasta_index

    script:
    """
    run_samtools.sh \
    -f ${fasta} \
    -n \
    ${fasta.name} \
    """

}
