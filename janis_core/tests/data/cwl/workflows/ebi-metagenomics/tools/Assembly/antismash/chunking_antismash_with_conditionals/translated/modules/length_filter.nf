nextflow.enable.dsl=2

process LENGTH_FILTER {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/before_qc/length_filter"
    cpus "${params.before_qc.length_filter.cpus}"
    memory "${params.before_qc.length_filter.memory}"

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
    ${params.before_qc.length_filter_stats_file_name} \
    ${submitted_seq_count} \
    --min_length ${min_length} \
    --extension ${params.before_qc.length_filter_input_file_format} \
    """

}
