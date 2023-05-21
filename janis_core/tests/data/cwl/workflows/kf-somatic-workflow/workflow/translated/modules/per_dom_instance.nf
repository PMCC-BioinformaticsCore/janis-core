nextflow.enable.dsl=2

process PER_DOM_INSTANCE {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/chop_and_avg_for_pfam2_cath/chop_and_avg_from_list/per_dom_instance"

    input:
    path fam_structs, stageAs: 'fam_structs'

    output:
    path "*.json", emit: dom_per_fam
    stdout, emit: family_name

    script:
    """
    python3 crossmapped_per_unp_dom.py \
    -f ${fam_structs} \
    """

}
