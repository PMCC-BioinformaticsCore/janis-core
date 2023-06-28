nextflow.enable.dsl=2

process ALIGN_AVG_STRUCTS_PAIRWISE {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/align_avg_structs_pairwise"

    input:
    val cath_fam_avg
    val core_avg
    val pfam_fam_avg

    output:
    path "<js>${ if (typeof inputs.result_file === 'string') {return inputs.result_file} else {return [ inputs.result_file.basename]}}</js>", emit: alignment_out

    script:
    """
    python3 pairwise_aligner.py \
    -c ${cath_fam_avg} \
    -p ${pfam_fam_avg} \
    -t ${core_avg} \
    -r align_Struct_analysis.csv \
    """

}
