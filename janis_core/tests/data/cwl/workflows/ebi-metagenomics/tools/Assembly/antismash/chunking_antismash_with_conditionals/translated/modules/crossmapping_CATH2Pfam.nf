nextflow.enable.dsl=2

process CROSSMAPPING_CATH2PFAM {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/crossmapping_cath2pfam"

    input:
    path cath_unq, stageAs: 'cath_unq'

    output:
    path "<js>${ if (typeof inputs.crossmap_cath === 'string') {return inputs.crossmap_cath} else {return [ inputs.crossmap_cath.basename]}}</js>", emit: allcrossmap_cath
    path "*.json", emit: cath_crossmapped
    path "<js>${ if (typeof inputs.no_crossmap === 'string') {return inputs.no_crossmap} else {return [ inputs.no_crossmap.basename]}}</js>", emit: cath_unmapped

    script:
    def min_dom_len = params.min_domain_length ? params.min_domain_length : 31
    """
    python3 map_unique_struct_cath2pfam.py \
    -l ${min_dom_len} \
    -p ../Data/pdbmap \
    -u cath_unq_unmapped.jsonx \
    -x cath_crossMapped_pfam.jsonx \
    -c ${cath_unq} \
    """

}
