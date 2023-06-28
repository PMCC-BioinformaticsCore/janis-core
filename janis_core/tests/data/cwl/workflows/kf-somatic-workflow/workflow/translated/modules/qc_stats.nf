nextflow.enable.dsl=2

process QC_STATS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/before_qc/qc_stats"

    input:
    path qced_reads, stageAs: 'qced_reads'
    val sequence_count

    output:
    path "inputs.out_dir_name", emit: output_dir
    path "<js>inputs.out_dir_name)/$(inputs.summary</js>", emit: summary_out

    script:
    """
    MGRAST_base.py \
    -i ${qced_reads} \
    -o <js>inputs.out_dir_name)/$(inputs.summary</js> \
    -d <js>${ var suffix = '.full';   if (inputs.sequence_count > inputs.max_seq) {     suffix = '.sub-set';   }   return "".concat(inputs.out_dir_name, '/', inputs.nucleotide_distribution, suffix);}</js> \
    -g <js>${ var suffix = '.full';   if (inputs.sequence_count > inputs.max_seq) {     suffix = '.sub-set';   }   return "".concat(inputs.out_dir_name, '/', inputs.gc_sum, suffix);}</js> \
    -l <js>${ var suffix = '.full';   if (inputs.sequence_count > inputs.max_seq) {     suffix = '.sub-set';   }   return "".concat(inputs.out_dir_name, '/', inputs.length_sum, suffix);}</js> \
    <js>${ if (inputs.sequence_count > inputs.max_seq) { return '-m '.concat(inputs.max_seq)} else { return ''} }</js> \
    """

}


process QC_STATS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/before_qc/qc_stats"

    input:
    path qced_reads, stageAs: 'qced_reads'
    val sequence_count

    output:
    path "inputs.out_dir_name", emit: output_dir
    path "<js>inputs.out_dir_name)/$(inputs.summary</js>", emit: summary_out

    script:
    """
    MGRAST_base.py \
    -i ${qced_reads} \
    -o <js>inputs.out_dir_name)/$(inputs.summary</js> \
    -d <js>${ var suffix = '.full';   if (inputs.sequence_count > inputs.max_seq) {     suffix = '.sub-set';   }   return "".concat(inputs.out_dir_name, '/', inputs.nucleotide_distribution, suffix);}</js> \
    -g <js>${ var suffix = '.full';   if (inputs.sequence_count > inputs.max_seq) {     suffix = '.sub-set';   }   return "".concat(inputs.out_dir_name, '/', inputs.gc_sum, suffix);}</js> \
    -l <js>${ var suffix = '.full';   if (inputs.sequence_count > inputs.max_seq) {     suffix = '.sub-set';   }   return "".concat(inputs.out_dir_name, '/', inputs.length_sum, suffix);}</js> \
    <js>${ if (inputs.sequence_count > inputs.max_seq) { return '-m '.concat(inputs.max_seq)} else { return ''} }</js> \
    """

}


process QC_STATS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/before_qc/qc_stats"

    input:
    path qced_reads, stageAs: 'qced_reads'
    val sequence_count

    output:
    path "inputs.out_dir_name", emit: output_dir
    path "<js>inputs.out_dir_name)/$(inputs.summary</js>", emit: summary_out

    script:
    """
    MGRAST_base.py \
    -i ${qced_reads} \
    -o <js>inputs.out_dir_name)/$(inputs.summary</js> \
    -d <js>${ var suffix = '.full';   if (inputs.sequence_count > inputs.max_seq) {     suffix = '.sub-set';   }   return "".concat(inputs.out_dir_name, '/', inputs.nucleotide_distribution, suffix);}</js> \
    -g <js>${ var suffix = '.full';   if (inputs.sequence_count > inputs.max_seq) {     suffix = '.sub-set';   }   return "".concat(inputs.out_dir_name, '/', inputs.gc_sum, suffix);}</js> \
    -l <js>${ var suffix = '.full';   if (inputs.sequence_count > inputs.max_seq) {     suffix = '.sub-set';   }   return "".concat(inputs.out_dir_name, '/', inputs.length_sum, suffix);}</js> \
    <js>${ if (inputs.sequence_count > inputs.max_seq) { return '-m '.concat(inputs.max_seq)} else { return ''} }</js> \
    """

}
