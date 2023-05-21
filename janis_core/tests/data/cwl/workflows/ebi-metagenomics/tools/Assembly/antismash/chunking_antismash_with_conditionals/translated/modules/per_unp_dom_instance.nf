nextflow.enable.dsl=2

process PER_UNP_DOM_INSTANCE {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/unmapped_from_pfam/per_unp_dom_instance"

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
