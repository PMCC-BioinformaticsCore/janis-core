nextflow.enable.dsl=2

process COMPARE_INSTANCES_CATH_PFAM {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/compare_instances_cath_pfam"

    input:
    path resmapped_cath, stageAs: 'resmapped_cath'
    path resmapped_pfam, stageAs: 'resmapped_pfam'
    path truedomains_file, stageAs: 'truedomains_file'

    output:
    path "inputs.unique_cath_struct", emit: cath_unique
    path "<js>${ if (typeof inputs.truedomains_file === 'string') {return inputs.truedomains_file} else {return [ inputs.truedomains_file.basename]}}</js>", emit: common_domains
    path "inputs.unique_pfam_struct", emit: pfam_unique

    script:
    def min_dom_len = params.min_domain_length ? params.min_domain_length : 31
    def truedomains_file = truedomains_file ? "-f ${truedomains_file}" : ""
    """
    python3 compare_cath_pfam.py \
    -l ${min_dom_len} \
    -c ${resmapped_cath} \
    -p ${resmapped_pfam} \
    ${truedomains_file} \
    -uq_pf unique_pfam.csv \
    -uq_ca unique_cath.csv \
    """

}
