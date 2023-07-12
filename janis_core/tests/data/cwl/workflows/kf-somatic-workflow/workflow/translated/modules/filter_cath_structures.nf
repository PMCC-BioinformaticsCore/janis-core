nextflow.enable.dsl=2

process FILTER_CATH_STRUCTURES {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/cath_domain_instances/filter_cath_structures"

    input:
    path track_fams, stageAs: 'track_fams'
    val min_dom_len

    output:
    path "<js>${ if (typeof inputs.obsolete_cath === 'string') {return inputs.obsolete_cath} else {return [ inputs.obsolete_cath.basename]}}</js>", emit: cath_obs
    path "<js>${ if (typeof inputs.separate_cath === 'string') {return inputs.separate_cath} else {return [ inputs.separate_cath.basename]}}</js>", emit: cath_structs
    path "*{inputs.split_suffix}", emit: splitted_cath_sep

    script:
    def min_dom_len = min_dom_len ? min_dom_len : 31
    """
    python3 separate_cath.py \
    -l ${min_dom_len} \
    -c ../Data/cath-domain-description-file.txt \
    -d ../Data/obsolete_PDB_entry_ids.txt \
    -n Filtered_CATH.csv \
    -o obsolete_cath.txt \
    -s part.csv \
    -f ${track_fams} \
    """

}
