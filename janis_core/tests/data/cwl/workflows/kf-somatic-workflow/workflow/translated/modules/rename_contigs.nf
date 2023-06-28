nextflow.enable.dsl=2

process RENAME_CONTIGS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/antismash/chunking/rename_contigs"

    input:
    path chunks, stageAs: 'chunks'
    path full_fasta, stageAs: 'full_fasta'
    val accession

    output:
    path "{inputs.chunks.basename}.tbl", emit: names_table
    path "antismash.*", emit: renamed_contigs_in_chunks

    script:
    """
    antismash_rename_contigs.py \
    -c ${chunks} \
    -i ${full_fasta} \
    -a ${accession} \
    """

}
