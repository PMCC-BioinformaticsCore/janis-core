nextflow.enable.dsl=2

process GEMSTATS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/gemstats"

    input:
    path carveme_gems

    output:
    path "{inputs.identifier}_CarveMe_GEMstats.tsv", emit: carveme_GEMstats

    script:
    def carveme_gems = carveme_gems.join(' ')
    """
    bash -x script.sh \
    ${params.identifier} \
    ${carveme_gems} \
    """

}
