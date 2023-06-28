nextflow.enable.dsl=2

process CHECK_ALIGNMENT_SCORES {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/check_alignment_scores"

    input:
    path aln_res_file, stageAs: 'aln_res_file'
    path fam_tracker, stageAs: 'fam_tracker'
    path cath_crossmap, stageAs: 'cath_crossmap'
    path pfam_crossmap, stageAs: 'pfam_crossmap'

    output:
    path "crossmapped_cath_failed.json", optional: true, emit: cath_crossmap_failed
    path "crossmapped_cath_passed.json", optional: true, emit: cath_crossmap_passed
    path "crossmapped_pfam_failed.json", optional: true, emit: pfam_crossmap_failed
    path "crossmapped_pfam_passed.json", optional: true, emit: pfam_crossmap_passed

    script:
    def aln_score = params.alignment_score ? params.alignment_score : Mscore
    def cath_crossmap = cath_crossmap ? "-cx ${cath_crossmap}" : ""
    def pfam_crossmap = pfam_crossmap ? "-px ${pfam_crossmap}" : ""
    def threshold_val = params.score_threshold ? params.score_threshold : 0.6
    """
    python3 check_threshold.py \
    -a ${aln_res_file} \
    -f ${fam_tracker} \
    -s ${aln_score} \
    -t ${threshold_val} \
    ${pfam_crossmap} \
    ${cath_crossmap} \
    """

}
