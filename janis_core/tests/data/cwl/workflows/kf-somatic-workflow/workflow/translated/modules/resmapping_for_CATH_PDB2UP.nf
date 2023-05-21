nextflow.enable.dsl=2

process RESMAPPING_FOR_CATH_PDB2_UP {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/cath_domain_instances/resmapping_cath_structs/resmapping_for_cath_pdb2_up"

    input:
    path cath_sep, stageAs: 'cath_sep'
    path sifts_dir, stageAs: 'sifts_dir'

    output:
    path "<js>${ if (typeof inputs.reslost === 'string') {return inputs.reslost} else {return [ inputs.reslost.basename]}}</js>", emit: cath_lost
    path "<js>${ if (typeof inputs.resmapping_file === 'string') {return inputs.resmapping_file} else {return [ inputs.resmapping_file.basename]}}</js>", emit: cath_resmapped

    script:
    """
    python3 resmapping_cath2up.py \
    -f ${cath_sep} \
    -s ${sifts_dir} \
    -m cath_resMapped.csv \
    -l lost_cath.txt \
    """

}
