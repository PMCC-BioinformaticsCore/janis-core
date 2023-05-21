nextflow.enable.dsl=2

process CHECK_THRESHOLD_STEP {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/unmapped_from_pfam/check_threshold_step"

    input:
    path all_struct_list, stageAs: 'all_struct_list'
    path aln_result, stageAs: 'aln_result'
    path failed_structs
    path passed_structs
    val aln_score
    val threshold_val

    output:
    path "<js>${ if (typeof inputs.failed_structs === 'string') {return inputs.failed_structs} else {return [ inputs.failed_structs.basename]}}</js>", emit: failed_structs_list
    path "<js>${ if (typeof inputs.passed_structs === 'string') {return inputs.passed_structs} else {return [ inputs.passed_structs.basename]}}</js>", emit: passed_structs_list

    script:
    def aln_score = aln_score ? aln_score : Mscore
    def failed_structs = failed_structs ? failed_structs : failed_structures_list.csv
    def passed_structs = passed_structs ? passed_structs : passed_structures_list.csv
    def threshold_val = threshold_val ? threshold_val : 0.6
    """
    python3 filter_align_scores.py \
    -f ${failed_structs} \
    -p ${passed_structs} \
    -s ${aln_score} \
    -t ${threshold_val} \
    -i ${aln_result} \
    -x ${all_struct_list} \
    """

}
