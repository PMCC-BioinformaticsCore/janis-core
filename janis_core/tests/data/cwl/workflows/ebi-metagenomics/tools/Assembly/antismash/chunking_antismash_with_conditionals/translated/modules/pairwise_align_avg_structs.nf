nextflow.enable.dsl=2

process PAIRWISE_ALIGN_AVG_STRUCTS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/unmapped_from_pfam/pairwise_align_avg_structs"

    input:
    val core_avg
    path query_dir, stageAs: 'query_dir'

    output:
    path "<js>${ if (typeof inputs.result_file === 'string') {return inputs.result_file} else {return [ inputs.result_file.basename]}}</js>", emit: alignment_out

    script:
    def query_dir = query_dir ? "-d ${query_dir}" : ""
    """
    python3 pairwise_aligner.py \
    ${query_dir} \
    -t ${core_avg} \
    -r align_Struct_analysis.csv \
    """

}
