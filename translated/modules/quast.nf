nextflow.enable.dsl=2

process QUAST {
    
    container "ppp-janis-translate:quast-5.0.2"
    publishDir "${params.outdir}/quast"

    input:
    path unknown1

    output:
    path "outputdir/report.html", emit: outReportHtml

    script:
    """
    quast metaquast \
    """

}
