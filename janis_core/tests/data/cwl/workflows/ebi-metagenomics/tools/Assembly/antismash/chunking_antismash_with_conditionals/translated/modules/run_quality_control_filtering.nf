nextflow.enable.dsl=2

process RUN_QUALITY_CONTROL_FILTERING {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/before_qc/run_quality_control_filtering"
    cpus "${params.before_qc.run_quality_control_filtering.cpus}"
    memory "${params.before_qc.run_quality_control_filtering.memory}"

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
    ${params.before_qc.run_quality_control_filtering_stats_file_name} \
    ${submitted_seq_count} \
    --min_length ${min_length} \
    --extension ${params.before_qc.run_quality_control_filtering_input_file_format} \
    """

}
