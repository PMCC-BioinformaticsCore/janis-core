nextflow.enable.dsl=2

process FILTER_PFAM_STRUCTURES {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/pfam_domain_instances/filter_pfam_structures"

    input:
    path track_fams, stageAs: 'track_fams'
    val min_dom_len

    output:
    path "<js>${ if (typeof inputs.obsolete_pfam === 'string') {return inputs.obsolete_pfam} else {return [ inputs.obsolete_pfam.basename]}}</js>", emit: pfam_obs
    path "<js>${ if (typeof inputs.separate_pfam === 'string') {return inputs.separate_pfam} else {return [ inputs.separate_pfam.basename]}}</js>", emit: pfam_structs
    path "*{inputs.split_suffix}", emit: splitted_pfam_sep

    script:
    def min_dom_len = min_dom_len ? min_dom_len : 31
    """
    python3 separate_pfam.py \
    -l ${min_dom_len} \
    -d ../Data/obsolete_PDB_entry_ids.txt \
    -p ../Data/pdbmap \
    -n Filtered_Pfam.csv \
    -o obsolete_pfam.txt \
    -s part.csv \
    -f ${track_fams} \
    """

}
