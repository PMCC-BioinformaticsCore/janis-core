nextflow.enable.dsl=2

process GET_FAMILY_IDS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/get_family_ids"

    input:
    path fam_tracker, stageAs: 'fam_tracker'

    output:
    path "<js>${ if (typeof inputs.fam_tracker === 'string') {return inputs.fam_tracker} else {return [ inputs.fam_tracker.basename]}}</js>", emit: family_ids

    script:
    def cath_ids = params.cath ? "-c " + params.cath.join(' ') : ""
    def fam_tracker = fam_tracker ? fam_tracker : family_ids.json
    def pfam_ids = params.pfam ? "-p " + params.pfam.join(' ') : ""
    """
    python3 get_family_ids.py \
    ${pfam_ids} \
    ${cath_ids} \
    -n ${params.iteration} \
    -f ${fam_tracker} \
    """

}
