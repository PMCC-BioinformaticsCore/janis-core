nextflow.enable.dsl=2

process ADD_CROSSMAPPED_TO_RESMAPPED {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/add_crossmapped_to_resmapped"

    input:
    path cath_crossmapped, stageAs: 'cath_crossmapped'
    path cath_resmapped, stageAs: 'cath_resmapped'
    path pfam_crossmapped, stageAs: 'pfam_crossmapped'
    path pfam_resmapped, stageAs: 'pfam_resmapped'

    output:
    path "<js>${ if (typeof inputs.cath_result === 'string') {return inputs.cath_result} else {return inputs.cath_result.basename}}</js>", emit: cath_structs
    path "<js>${ if (typeof inputs.pfam_result === 'string') {return inputs.pfam_result} else {return inputs.pfam_result.basename}}</js>", emit: pfam_structs

    script:
    def cath_crossmapped = cath_crossmapped ? "-cx ${cath_crossmapped}" : ""
    def cath_resmapped = cath_resmapped ? "-c ${cath_resmapped}" : ""
    def pfam_crossmapped = pfam_crossmapped ? "-px ${pfam_crossmapped}" : ""
    def pfam_resmapped = pfam_resmapped ? "-p ${pfam_resmapped}" : ""
    """
    python3 add_crossmapped2resmapped.py \
    ${pfam_resmapped} \
    ${cath_resmapped} \
    ${pfam_crossmapped} \
    ${cath_crossmapped} \
    -pr pfam_res_crossMapped.csv \
    -cr cath_res_crossMapped.csv \
    """

}
