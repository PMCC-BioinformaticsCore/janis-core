nextflow.enable.dsl=2

process BUSCO {
    
    container "ppp-janis-translate:busco-5.3.2"
    publishDir "${params.outdir}/busco"

    input:
    path in_file

    output:
    path "busco_galaxy/run_*/short_summary.txt", emit: outBuscoSum
    path "busco_galaxy/run_*/full_table.tsv", emit: outBuscoTable

    script:
    """
    busco \
    --in ${in_file} \
    """

}
