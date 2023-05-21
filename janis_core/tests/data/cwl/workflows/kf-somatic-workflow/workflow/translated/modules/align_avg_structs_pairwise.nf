nextflow.enable.dsl=2

process ALIGN_AVG_STRUCTS_PAIRWISE {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/align_avg_structs_pairwise"

    input:
    path core_avg
    path cath_fam_avg
    path pfam_fam_avg

    output:
    path "<js>${ if (typeof inputs.result_file === 'string') {return inputs.result_file} else {return [ inputs.result_file.basename]}}</js>", emit: alignment_out

    script:
    def cath_fam_avg = cath_fam_avg ? "-c " + cath_fam_avg.join(' ') : ""
    def core_avg = core_avg.join(' ')
    def pfam_fam_avg = pfam_fam_avg ? "-p " + pfam_fam_avg.join(' ') : ""
    """
    python3 pairwise_aligner.py \
    ${cath_fam_avg} \
    ${pfam_fam_avg} \
    -t ${core_avg} \
    -r align_Struct_analysis.csv \
    """

}
