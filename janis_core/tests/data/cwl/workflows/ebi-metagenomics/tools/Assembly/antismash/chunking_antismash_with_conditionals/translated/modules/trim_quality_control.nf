nextflow.enable.dsl=2

process TRIM_QUALITY_CONTROL {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/before_qc/trim_quality_control"
    cpus "${params.before_qc.trim_quality_control.cpus}"
    memory "${params.before_qc.trim_quality_control.memory}"

    input:
    path reads1, stageAs: 'reads1'

    output:
    path "trim.log", emit: log_file
    path "{inputs.reads1.nameroot}.trimmed", emit: reads1_trimmed
    path "{inputs.reads1.basename}.trimmed.unpaired.fastq", optional: true, emit: reads1_trimmed_unpaired

    script:
    def leading = params.before_qc.trim_quality_control_leading ? "LEADING:${params.before_qc.trim_quality_control_leading}" : ""
    def minlen = params.before_qc.trim_quality_control_minlen ? "MINLEN:${params.before_qc.trim_quality_control_minlen}" : ""
    def phred = params.before_qc.trim_quality_control_phred ? "-phred${params.before_qc.trim_quality_control_phred}" : ""
    def slidingwindow = params.before_qc.trim_quality_control_slidingwindow ? "SLIDINGWINDOW:${params.before_qc.trim_quality_control_slidingwindow}" : ""
    def trailing = params.before_qc.trim_quality_control_trailing ? "TRAILING:${params.before_qc.trim_quality_control_trailing}" : ""
    """
    trimmomatic.sh \
    ${params.before_qc.trim_quality_control_end_mode} \
    ${phred} \
    -threads <js>runtime.cores</js> \
    -trimlog trim.log \
    ${reads1} \
    "${${reads1}.baseName}.trimmed" \
    ${leading} \
    ${trailing} \
    ${slidingwindow} \
    ${minlen} \
    """

}
