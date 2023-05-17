nextflow.enable.dsl=2

process AVG_AVERAGED_STRUCTURES {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/chop_and_avg_for_core/avg_averaged_structures"

    input:
    path fam_name, stageAs: 'fam_name'
    path split_dir, stageAs: 'split_dir'

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
