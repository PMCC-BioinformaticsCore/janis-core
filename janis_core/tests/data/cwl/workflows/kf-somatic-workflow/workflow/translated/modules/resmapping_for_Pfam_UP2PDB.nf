nextflow.enable.dsl=2

process RESMAPPING_FOR_PFAM_UP2_PDB {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/pfam_domain_instances/resmapping_pfam_structs/resmapping_for_pfam_up2_pdb"

    input:
    path pfam_sep, stageAs: 'pfam_sep'
    path sifts_dir, stageAs: 'sifts_dir'

    output:
    path "<js>${ if (typeof inputs.reslost === 'string') {return inputs.reslost} else {return [ inputs.reslost.basename]}}</js>", emit: pfam_lost
    path "<js>${ if (typeof inputs.resmapping_file === 'string') {return inputs.resmapping_file} else {return [ inputs.resmapping_file.basename]}}</js>", emit: pfam_resmapped

    script:
    """
    python3 resmapping_pfam2pdb.py \
    -f ${pfam_sep} \
    -s ${sifts_dir} \
    -m pfam_resMapped.csv \
    -l lost_pfam.txt \
    """

}
