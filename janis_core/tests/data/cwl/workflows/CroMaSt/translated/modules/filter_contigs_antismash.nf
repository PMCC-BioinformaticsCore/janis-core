nextflow.enable.dsl=2

process FILTER_CONTIGS_ANTISMASH {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/antismash/filtering/filter_contigs_antismash"
    cpus "${params.after_qc.antismash.filtering.filter_contigs_antismash.cpus}"
    memory "${params.after_qc.antismash.filtering.filter_contigs_antismash.memory}"

    input:
    path seq_file, stageAs: 'seq_file'
    val min_length
    val submitted_seq_count

    output:
    path "{inputs.seq_file.nameroot}.fasta", emit: filtered_file
    path "inputs.stats_file_name", emit: stats_summary_file

    script:
    """
    run_quality_filtering.py \
    ${seq_file} \
    "${${seq_file}.baseName}.fasta" \
    ${params.after_qc.antismash.filtering.filter_contigs_antismash_stats_file_name} \
    ${submitted_seq_count} \
    --min_length ${min_length} \
    --extension ${params.after_qc.antismash.filtering.filter_contigs_antismash_input_file_format} \
    """

}
