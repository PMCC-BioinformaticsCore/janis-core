nextflow.enable.dsl=2

process GATHER_FAILED_DOMAINS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/gather_failed_domains"

    input:
    path cath_crossmapped, stageAs: 'cath_crossmapped'
    path cath_unmapped, stageAs: 'cath_unmapped'
    path pfam_crossmapped, stageAs: 'pfam_crossmapped'
    path pfam_unmapped, stageAs: 'pfam_unmapped'
    path unmapped_out, stageAs: 'unmapped_out'

    output:
    path "<js>${ if (typeof inputs.unmapped_out === 'string') {return inputs.unmapped_out} else {return [ inputs.unmapped_out.basename]}}</js>", emit: unmapped_list

    script:
    def cath_crossmapped = cath_crossmapped ? "-cx ${cath_crossmapped}" : ""
    def cath_unmapped = cath_unmapped ? "-c ${cath_unmapped}" : ""
    def pfam_crossmapped = pfam_crossmapped ? "-px ${pfam_crossmapped}" : ""
    def pfam_unmapped = pfam_unmapped ? "-p ${pfam_unmapped}" : ""
    def unmapped_out = unmapped_out ? unmapped_out : domain_StIs_f.json
    """
    python3 merge_unmapped.py \
    ${pfam_unmapped} \
    ${cath_unmapped} \
    -o ${unmapped_out} \
    ${pfam_crossmapped} \
    ${cath_crossmapped} \
    """

}
