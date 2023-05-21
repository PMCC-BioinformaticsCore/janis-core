nextflow.enable.dsl=2

process ADD_DOMAIN_POSITIONS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/pfam_domain_instances/add_domain_positions"

    input:
    path resmapped_files
    path dom_posi_file, stageAs: 'dom_posi_file'

    output:
    path "<js>${ if (typeof inputs.dom_posi_file === 'string') {return inputs.dom_posi_file} else {return [ inputs.dom_posi_file.basename ]}}</js>", emit: resmapped_domains

    script:
    def dom_posi_file = dom_posi_file ? dom_posi_file : resmapped_domains_posi.csv
    def resmapped_files = resmapped_files.join(' ')
    """
    python3 add_domain_num.py \
    -i ${resmapped_files} \
    -o ${dom_posi_file} \
    """

}
