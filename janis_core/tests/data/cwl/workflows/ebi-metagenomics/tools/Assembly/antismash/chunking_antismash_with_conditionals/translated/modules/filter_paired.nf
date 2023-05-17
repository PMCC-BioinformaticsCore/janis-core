nextflow.enable.dsl=2

process FILTER_PAIRED {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/before_qc/overlap_reads/filter_paired"
    cpus "${params.before_qc.overlap_reads.filter_paired.cpus}"
    memory "${params.before_qc.overlap_reads.filter_paired.memory}"

    input:
    path fastq1, stageAs: 'fastq1'
    path fastq2, stageAs: 'fastq2'
    val min_length_required

    output:
    path "fastp.html", emit: html_report
    path "fastp.json", emit: json_report
    path "{inputs.fastq1.nameroot}.fastp.fastq", emit: out_fastq1
    path "{inputs.fastq2.nameroot}.fastp.fastq", optional: true, emit: out_fastq2

    script:
    def base_correction = params.before_qc.overlap_reads.filter_paired_base_correction == false ? "" : "--correction"
    def disable_trim_poly_g = params.before_qc.overlap_reads.filter_paired_disable_trim_poly_g == false ? "" : "--disable_trim_poly_g"
    def fastq2 = fastq2 ? "-I ${fastq2}" : ""
    def force_polyg_tail_trimming = params.before_qc.overlap_reads.filter_paired_force_polyg_tail_trimming == false ? "" : "--trim_poly_g"
    def min_length_required = min_length_required ? min_length_required : 50
    def threads = params.before_qc.overlap_reads.filter_paired_threads ? params.before_qc.overlap_reads.filter_paired_threads : 1
    """
    fastp \
    ${fastq2} \
    -i ${fastq1} \
    --length_required ${min_length_required} \
    --thread ${threads} \
    ${base_correction} \
    ${disable_trim_poly_g} \
    ${force_polyg_tail_trimming} \
    --qualified_quality_phred 20 \
    --unqualified_percent_limit 20 \
    -o "${${fastq1}.baseName}.fastp.fastq" \
    <js>${  if (inputs.fastq2){    return '-O';  } else {    return '';  }}</js> \
    <js>${  if (inputs.fastq2){    return inputs.fastq2.nameroot + ".fastp.fastq";  } else {    return '';  }}</js> \
    """

}


process FILTER_PAIRED {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/before_qc/overlap_reads/filter_paired"
    cpus "${params.before_qc.overlap_reads.filter_paired.cpus}"
    memory "${params.before_qc.overlap_reads.filter_paired.memory}"

    input:
    path fastq1, stageAs: 'fastq1'
    path fastq2, stageAs: 'fastq2'
    val min_length_required

    output:
    path "fastp.html", emit: html_report
    path "fastp.json", emit: json_report
    path "{inputs.fastq1.nameroot}.fastp.fastq", emit: out_fastq1
    path "{inputs.fastq2.nameroot}.fastp.fastq", optional: true, emit: out_fastq2

    script:
    def base_correction = params.before_qc.overlap_reads.filter_paired_base_correction == false ? "" : "--correction"
    def disable_trim_poly_g = params.before_qc.overlap_reads.filter_paired_disable_trim_poly_g == false ? "" : "--disable_trim_poly_g"
    def fastq2 = fastq2 ? "-I ${fastq2}" : ""
    def force_polyg_tail_trimming = params.before_qc.overlap_reads.filter_paired_force_polyg_tail_trimming == false ? "" : "--trim_poly_g"
    def min_length_required = min_length_required ? min_length_required : 50
    def threads = params.before_qc.overlap_reads.filter_paired_threads ? params.before_qc.overlap_reads.filter_paired_threads : 1
    """
    fastp \
    ${fastq2} \
    -i ${fastq1} \
    --length_required ${min_length_required} \
    --thread ${threads} \
    ${base_correction} \
    ${disable_trim_poly_g} \
    ${force_polyg_tail_trimming} \
    --qualified_quality_phred 20 \
    --unqualified_percent_limit 20 \
    -o "${${fastq1}.baseName}.fastp.fastq" \
    <js>${  if (inputs.fastq2){    return '-O';  } else {    return '';  }}</js> \
    <js>${  if (inputs.fastq2){    return inputs.fastq2.nameroot + ".fastp.fastq";  } else {    return '';  }}</js> \
    """

}
