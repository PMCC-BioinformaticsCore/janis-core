nextflow.enable.dsl=2

process COLLECT_LOST_INSTANCES {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/pfam_domain_instances/collect_lost_instances"

    input:
    path lost_instance
    path obs_insta, stageAs: 'obs_insta'
    path outfile

    output:
    path "<js>${ if (typeof inputs.outfile === 'string') {return inputs.outfile} else {return inputs.outfile.basename}}</js>", emit: lost_domain_list

    script:
    def outfile = outfile ? outfile : lost_resmap.json
    """
    python3 collect_lost_instances.py \
    ${outfile} \
    """

}
