nextflow.enable.dsl=2

process FORMAT_CORE_LIST {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/format_core_list"

    input:
    path infile, stageAs: 'infile'

    output:
    path "<js>${ if (typeof inputs.outfile === 'string') {return inputs.outfile} else {return inputs.outfile.basename}}</js>", emit: coredomains_list

    script:
    """
    python3 list_true_domains.py \
    """

}
