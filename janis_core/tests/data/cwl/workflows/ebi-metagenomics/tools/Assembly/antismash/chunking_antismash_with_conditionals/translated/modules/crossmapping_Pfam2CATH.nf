nextflow.enable.dsl=2

process CROSSMAPPING_PFAM2_CATH {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/crossmapping_pfam2_cath"

    input:
    path pfam_unq, stageAs: 'pfam_unq'

    output:
    path "<js>${ if (typeof inputs.crossmap_pfam === 'string') {return inputs.crossmap_pfam} else {return [ inputs.crossmap_pfam.basename]}}</js>", emit: allcrossmap_pfam
    path "*.json", emit: pfam_crossmapped
    path "<js>${ if (typeof inputs.no_crossmap === 'string') {return inputs.no_crossmap} else {return [ inputs.no_crossmap.basename]}}</js>", emit: pfam_unmapped

    script:
    def min_dom_len = params.min_domain_length ? params.min_domain_length : 31
    """
    python3 map_unique_struct_pfam2cath.py \
    -l ${min_dom_len} \
    -c ../Data/cath-domain-description-file.txt \
    -u pfam_unq_unmapped.jsonx \
    -x pfam_crossMapped_cath.jsonx \
    -p ${pfam_unq} \
    """

}
