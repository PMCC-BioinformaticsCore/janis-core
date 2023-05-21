nextflow.enable.dsl=2

process AVG_CHOPPED_STRUCTS_UNP_DOMAINS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/unmapped_from_pfam/avg_unp_domains/avg_chopped_structs_unp_domains"

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
