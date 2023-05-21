nextflow.enable.dsl=2

process AVG_CHOPPED_STRUCTURES_UNP_DOMAINS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/chop_and_avg_for_pfam2_cath/chop_and_avg_from_list/avg_unp_domains/avg_chopped_structures_unp_domains"

    input:
    path fam_name, stageAs: 'fam_name'
    path split_dir

    output:
    path "*.pdb", emit: avg_structs

    script:
    def split_dir = split_dir ? split_dir : split_PDB
    """
    python3 align_compute_avg.py \
    -f ${fam_name} \
    -s ${split_dir} \
    -k KPAX_RESULTS \
    """

}
