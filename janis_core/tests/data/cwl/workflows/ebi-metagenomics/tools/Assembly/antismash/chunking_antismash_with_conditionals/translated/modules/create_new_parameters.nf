nextflow.enable.dsl=2

process CREATE_NEW_PARAMETERS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/create_new_parameters"

    input:
    path core_domain_struct
    path crossmap_cath, stageAs: 'crossmap_cath'
    path crossmap_pfam, stageAs: 'crossmap_pfam'
    path domain_like, stageAs: 'domain_like'
    path failed_domains, stageAs: 'failed_domains'
    path fam_tracker, stageAs: 'fam_tracker'
    path in_paramfile, stageAs: 'in_paramfile'
    path true_domains, stageAs: 'true_domains'

    output:
    path "<js>${ if (typeof inputs.next_paramfile === 'string') {return inputs.next_paramfile} else {return inputs.next_paramfile.basename}}</js>", emit: next_parmfile

    script:
    """
    python3 create_param.py \
    -i ${in_paramfile} \
    -o new_param.yml \
    -px ${crossmap_pfam} \
    -cx ${crossmap_cath} \
    """

}
