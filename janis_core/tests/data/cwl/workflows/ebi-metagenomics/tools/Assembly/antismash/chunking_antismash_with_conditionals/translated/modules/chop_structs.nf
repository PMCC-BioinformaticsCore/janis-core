nextflow.enable.dsl=2

process CHOP_STRUCTS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/unmapped_from_pfam/avg_unp_domains/chop_structs"

    input:
    path struct_insta, stageAs: 'struct_insta'
    path pdb_dir, stageAs: 'pdb_dir'

    output:
    stdout, emit: family_name
    path "inputs.split_dir", emit: split_structs_dir

    script:
    def pdb_dir = pdb_dir ? pdb_dir : ../PDB_files/
    """
    python3 chop_struct2domains.py \
    -f ${struct_insta} \
    -p ${pdb_dir} \
    -s split_PDB \
    -k KPAX_RESULTS \
    """

}
